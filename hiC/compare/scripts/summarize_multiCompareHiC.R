# 4.0.2
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(multiHiCcompare)
library(BiocParallel)
library(edgeR)

args = commandArgs(trailingOnly=TRUE)
dir=args[1]

# dir="/mnt/isilon/sfgi/suc1/analyses/wells/hiC/compare/naiveT_stim/data/4000"

data_files=dir(dir, pattern="HiCcompare\\.data\\.Rdata", recursive=T)

############################ read loop file ###########################
loops = read_delim(
        dir(dir,pattern="mustache_ICE_loops\\.merged", full.names=T), delim="\t",
        col_names=c("chr_a","start_a","end_a","chr_b","start_b","end_b")
)

resol = loops %>% mutate(resol=end_a-start_a) %>% distinct(resol) %>% pull(resol)

loops = loops %>% 
select(chr=chr_a, region1=start_a, region2=start_b) %>% 
mutate(chr=gsub("chr","",chr)) %>% 
mutate(chr=case_when(
        chr=="X" ~ "23", 
        chr=="Y" ~ "24",
        T ~ chr
)) %>% 
mutate(chr=as.numeric(chr))

######################### function read_HiCcompare ############################
# data_file=data_files[1]
read_HiCcompare <- function(data_file, padj_cutoff=0.05, abslogFC_cutoff=1){
        samples=dir(dirname(file.path(dir,data_file)), pattern="*.txt")
        samples=sapply(samples, function(l){strsplit(l,"\\.")[[1]][1]}, USE.NAMES=F)

        load(file.path(dir,data_file))
        colnames(normIF) = c("chr","region1","region2","D", samples)

        global_sigDiff_results = lapply(
                results,
                function(df){
                        df %>% as_tibble() %>% 
                        filter(p.adj < padj_cutoff, abs(logFC) > abslogFC_cutoff) %>% 
                        left_join(
                                normIF %>% as_tibble() %>% 
                                select(-D)
                        )
                }
        )


        loop_normIF = normIF %>% 
        as_tibble() %>% 
        semi_join(loops)

        list(
                global_sigDiff_results=global_sigDiff_results,
                loop_normIF=loop_normIF
        )
}

######################### read_HiCcompare run #################################
HiCcompare_results = lapply(
        data_files,
        read_HiCcompare
)

global_sigDiff_results = lapply(
        HiCcompare_results,
        function(x){
                x[["global_sigDiff_results"]]
        }
)

contrast_names=c("Alpha_vs_Acinar", "Beta_vs_Acinar", "Beta_vs_Alpha") # change, make sure same as multiCompareHiC_byChr.R

global_sigDiff_results = lapply(
        contrast_names,
        function(contrast_name){
                do.call(
                        "bind_rows",
                        lapply(
                                global_sigDiff_results,
                                function(x){
                                        x[[contrast_name]]
                                }
                        )
                )
        }
)
names(global_sigDiff_results) = contrast_names

loop_normIF = do.call(
        "bind_rows",
        lapply(
                HiCcompare_results,
                function(x){
                        x[["loop_normIF"]]
                }
        )
)

save(global_sigDiff_results, loop_normIF, loops, resol, file=file.path(dir,"HiCcompare_sigDiff.Rdata"))

############################# loop differential ################################

#### helper functions from diff_loops
.get_dist_tables <- function(chr_table, max.pool) {
  # check for correct max.pool input
  if (!is.numeric(max.pool)) {
    stop("max.pool must be a numeric value between 0 and 1")
  }
  if (max.pool < 0 || max.pool > 1) {
    stop("max.pool must be between 0 and 1")
  }
  if (max.pool < 0.5) {
    warning("Setting max.pool < 0.5 may cause issues")
  }
  if (max.pool > 0.7) {
    warning("Setting max.pool > 0.7 may cause issues")
  }
  D <- sort(unique(chr_table$D)) 
  # use triangular number series to solve for number of pools
  x <- length(D)
  n <- (sqrt((8 * x) + 1) - 1) / 2
  n <- ceiling(n)
  pools <- rep(D[1:n], 1:n)[1:x]
  # add the pools column to the data.table corresponding to the distances
  id.df <- as.data.frame(cbind(D, pools))
  # get maximum distance
  dist_max <- ceiling(max.pool * nrow(id.df))
  # combine pools for everything at or above maximum distance
  pool_max <- id.df[dist_max,2] 
  id.df[id.df$pools >= pool_max,2] <- pool_max
  table_by_dist <- data.table::as.data.table(left_join(chr_table, id.df, by = c("D" = "D")))
  # split up tables by pool
  table_by_dist <- split(table_by_dist, table_by_dist$pools)
  # drop pools column
  table_by_dist <- lapply(table_by_dist, function(x) x[, pools := NULL])
  # check number of rows per table
  table_by_dist <- .check_tables(table_by_dist)
  return(table_by_dist)
}

.check_tables <- function(table_by_dist) {
  min.row <- 500
  # get row numbers
  row.nums <- lapply(table_by_dist, nrow)
  # find which tables have number of rows < min.row
  small_D <- which(row.nums < min.row)
  # check if any tables have less than min.row
  if (length(small_D) > 0) {
    while (length(small_D) > 0) {
      i <- small_D[1]
      # check if first distance has less than min.row
      if(i==length(table_by_dist)){ # if it reach to the end, forget it
        return(table_by_dist)
      }
      if (table_by_dist[[i]]$D[1] == 0 | i==1) { # if i==1, it is the first one
        new_table <- rbind(table_by_dist[[i]], table_by_dist[[i+1]])
        table_by_dist[[i]] <- new_table
        table_by_dist[[i+1]] <- NA
      } else{
           # otherwise combine with previous table
           new_table <- rbind(table_by_dist[[i-1]], table_by_dist[[i]])
           table_by_dist[[i-1]] <- new_table
           table_by_dist[[i]] <- NA
      }
      # remove NAs
      table_by_dist <- table_by_dist[!is.na(table_by_dist)]
      # update row numbers
      row.nums <- lapply(table_by_dist, nrow)
      # update which tables have number of rows < min.row
      small_D <- which(row.nums < min.row)
    }
    return(table_by_dist)
  } else {
    return(table_by_dist)
  }
}

.hictable2DGEList <- function(hic_table, covariates) {
  # convert IFs into a matrix
  IFs <- as.matrix(hic_table[, -c("chr", "region1", "region2", "D"), with = FALSE])
  # create DGEList 
  dge <- DGEList(counts = IFs, samples = covariates)
  return(dge)
}

.glm_reformat <- function(result, hic_table, p.method) {
  # create table of location info and p-value results
  result <- cbind(hic_table, result$table)
  colnames(result)[ncol(result)] <-"p.value"
  # adjust p-values
  result$p.adj <- p.adjust(result$p.value, method = p.method)
  return(result)
}

#### function diff_loops
diff_loops <- function (loop_normIF, design, contrast = NA, coef = NA, method = "QLFTest", M = 1, p.method = "fdr", parallel = FALSE, max.pool = 0.7, groups) {
    method <- match.arg(method, c("LRTest", "QLFTest", "Treat"), 
        several.ok = FALSE)
    # if (!normalized(hicexp2)) {
    #     warning("You should normalize the data before entering it into hic_glm")
    # }
    if ((is.na(contrast[1]) & is.na(coef[1])) || (!is.na(contrast[1]) & 
        !is.na(coef[1]))) {
        stop("You must enter a value for contrast or a coef, but not both")
    }
    # table_list <- split(hic_table(hicexp2), hic_table(hicexp2)$chr)
    # table_list <- lapply(table_list, .get_dist_tables, max.pool = max.pool)
    table_list <- .get_dist_tables(loop_normIF, max.pool = max.pool)
    # table_list <- do.call(c, table_list)
    group = data.frame(groups)
    colnames(group)="group"
    dge_list <- lapply(table_list, .hictable2DGEList, covariates = group)
    dge_list <- lapply(dge_list, edgeR::estimateDisp, design = design)
    fit <- lapply(dge_list, edgeR::glmQLFit, design = design)
    if (method == "QLFTest") {
        if (is.na(coef)) {
            result <- lapply(fit, edgeR::glmQLFTest, contrast = contrast)
        }
        else {
            result <- lapply(fit, edgeR::glmQLFTest, coef = coef)
        }
    }
    if (method == "LRTest") {
        if (is.na(coef)) {
            result <- lapply(fit, edgeR::glmLRT, contrast = contrast)
        }
        else {
            result <- lapply(parallel, fit, edgeR::glmLRT, coef = coef)
        }
    }
    if (method == "Treat") {
        if (is.na(coef)) {
            result <- lapply(fit, edgeR::glmTreat, contrast = contrast, lfc = M)
        }
        else {
            result <- lapply(fit, edgeR::glmTreat, coef = coef, lfc = M)
        }
    }

    result <- mapply(.glm_reformat, result, table_list, 
            MoreArgs = list(p.method = p.method), SIMPLIFY = FALSE)
    result <- data.table::rbindlist(result)
    result <- result %>% arrange(p.adj)
    result
}

#### run loop differential
# sample, groups (double check)
samples = colnames(loop_normIF %>% select(-chr, -region1, -region2, -D))
groups = sapply(samples, function(x){strsplit(x,"_")[[1]][1]})

# design
design = model.matrix(~0+groups)
rownames(design) = samples
colnames(design) = gsub("groups","",colnames(design))

# contrast (change) - make sure same as previous contrast_names
my.contrasts = makeContrasts(
        Alpha_vs_Acinar=Alpha - Acinar,
        Beta_vs_Acinar=Beta-Acinar,
        Beta_vs_Alpha=Beta-Alpha,
        levels=design
)


loop_sigDiff_results = lapply(
        colnames(my.contrasts),
        function(contrast_name){
                diff_loops(
                        loop_normIF, 
                        design=design, 
                        contrast=my.contrasts[,contrast_name], 
                        groups=groups
                ) %>% as_tibble()
        }
)

names(loop_sigDiff_results) = colnames(my.contrasts)

save(loop_sigDiff_results, global_sigDiff_results, loop_normIF, loops, resol, file=file.path(dir,"HiCcompare_sigDiff.Rdata"))

######################## check how many loop differential on global way and loop way ###############
lapply(
        global_sigDiff_results,
        function(df){
                semi_join(
                        df,
                        loops
                )
        }
)

lapply(
        loop_sigDiff_results,
        function(df){
                df %>% 
                filter(abs(logFC) >= 1, p.adj < 0.05)
        }
)


########################## write result ####################

write_diffLoop <- function(comp, DF_list, prefix){
        df = DF_list[[comp]] %>% 
        filter(abs(logFC) >= 1, p.adj < 0.05)
        df = df %>% 
        dplyr::rename(
                chr_a=chr,
                start_a=region1,
                start_b=region2
        ) %>% 
        mutate(
                chr_b=chr_a,
                end_a=start_a+resol,
                end_b=start_b+resol
        ) %>% 
        select(chr_a, start_a, end_a, chr_b, start_b, end_b,logFC, p.value, p.adj, !!samples) %>% 
        mutate_at(c("start_a","end_a","start_b","end_b"), function(x){format(x, scientific=F)}) %>% 
        mutate_at(c("chr_a","chr_b"), function(x){gsub("^","chr",x)})
        write.table(
                df, file=file.path(dir,paste0(prefix,".",comp,".txt")),
                sep="\t", row.names=F, quote=F
        )
}

lapply(
        names(loop_sigDiff_results),
        write_diffLoop,
        DF_list = loop_sigDiff_results,
        prefix="loop_sigDiff"
)

lapply(
        names(loop_sigDiff_results),
        write_diffLoop,
        DF_list = lapply(
                global_sigDiff_results,
                function(df){
                        semi_join(
                                df,
                                loops
                        )
                }
        ),
        prefix="globalLoop_sigDiff"
)

lapply(
        names(loop_sigDiff_results),
        write_diffLoop,
        DF_list = global_sigDiff_results,
        prefix="global_sigDiff"
)



