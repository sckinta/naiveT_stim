# 4.0.2
# https://bioconductor.org/packages/release/bioc/vignettes/multiHiCcompare/inst/doc/multiHiCcompare.html
# BiocManager::install("multiHiCcompare")
# BiocManager::install("S4Vectors")
# BiocManager::install("BiocParallel")
library(dplyr)
library(tidyr)
library(ggplot2)
library(multiHiCcompare)
library(BiocParallel)


args = commandArgs(trailingOnly=TRUE)
dir=args[1]
out_prefix=args[2]
numCores=args[3]

# dir="/mnt/isilon/sfgi/suc1/analyses/grant/hiC/compare/pancreaticCells/data/1M"
# out_prefix="/mnt/isilon/sfgi/suc1/analyses/grant/hiC/compare/pancreaticCells/data/1M/multiCompareHiC.1M"
# numCores=1

if(numCores!=1){
        register(SnowParam(workers = numCores), default = TRUE)
        paral=T
}else{
        paral=F
}


files=dir(dir, pattern="*.txt", full.names=T)

data = lapply(
        files,
        function(file){
                df <- read.table(file, header = FALSE)
                # The chromosome number should be entered as just the number
                # Chromosomes such as X, Y, etc. should be entered as 23, 24, etc
                df$V1 = gsub("chr","",df$V1)
                df$V4 = gsub("chr","",df$V4)
                df = df %>% as_tibble() %>% 
                filter(V1!="M", V4!="M") %>% 
                mutate(V1=case_when(
                        V1=="X" ~ "23", 
                        V1=="Y" ~ "24",
                        T ~ V1
                )) %>% 
                mutate(V4=case_when(
                        V4=="X" ~ "23", 
                        V4=="Y" ~ "24",
                        T ~ V4
                ))
                
                sparse <- HiCcompare::cooler2sparse(df) # it will report error if there is only one chr
                sparse # region1, region2 are start point (0-based)
        }
)

names(data) = sapply(
        basename(files),
        function(s){strsplit(s,"\\.")[[1]][1]}
)


data_chrs=lapply(
        1:24,
        function(chr){
                lapply(
                        data,
                        function(df){
                                df[[chr]]$chr=chr
                                df[[chr]][,c("chr","region1","region2","IF")]
                        }
                )
        }
)

groups = as.integer(
        factor(sapply(
        names(data),
        function(s){strsplit(s,"_")[[1]][1]}
))
)
groups = groups-1

normTest_by_chr <- function(chr){
        # make hicexp
        hicexp1 <- make_hicexp(data_list=data_chrs[[chr]], 
                               groups = groups, 
                               zero.p = 0.8, A.min = 5, filter = TRUE,
                               remove.regions = hg19_cyto)

# The zero.p option allows for filtering by the proportion of zero IFs for an interaction. The A.min allows for filtering by a minimum average IF value. These options can be used together or individually to filter your data. Filtering is important to remove interactions with lots of 0 IFs and low average expression.

        # normalize
        hicexp1 <- fastlo(hicexp1, verbose = FALSE, parallel = paral)
        # hicexp1 <- cyclic_loess(hicexp1, verbose = FALSE, parallel = FALSE)

        # filename=paste0("MD_",chr,".pdf")
        # pdf(file.path(dir,"plots", filename))
        # MD_hicexp(hicexp1)
        # dev.off()

        # hicexp1 <- hic_exactTest(hicexp1, p.method = 'fdr', parallel = FALSE)

        # hic_table(hicexp1) # normalized IF
        # results(hicexp1)
        # meta(hicexp1)

        d <- model.matrix(~factor(meta(hicexp1)$group))

        hicexp1 <- hic_glm(
                hicexp1, design = d, coef = 2, 
                method = "QLFTest", p.method = "fdr", parallel = paral
        )

        list(
                results=results(hicexp1),
                normIF=hic_table(hicexp1)
        )
}

results = lapply(
        1:23,
        normTest_by_chr
)

# test_results
test_results = do.call(
        "bind_rows",
        lapply(
                results,
                function(l){l[["results"]]}
        )
)

test_results$p.adj=p.adjust(test_results$p.value, method = "fdr")


# normIF
normIF = do.call(
        "bind_rows",
        lapply(
                results,
                function(l){l[["normIF"]]}
        )
)

save(normIF, test_results, file=paste0(out_prefix,".data.Rdata"))

write.csv(
        test_results %>% as_tibble() %>% filter(abs(logFC) >1, p.adj < 0.05),
        file=paste0(out_prefix,".sigResults.csv"),
        row.names=F
)