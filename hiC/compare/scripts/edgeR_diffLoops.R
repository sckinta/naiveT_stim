function (hicexp2, design, contrast = NA, coef = NA, method = "QLFTest", 
    M = 1, p.method = "fdr", parallel = FALSE, max.pool = 0.7) 
{
    method <- match.arg(method, c("LRTest", "QLFTest", "Treat"), 
        several.ok = FALSE)
    if (!normalized(hicexp2)) {
        warning("You should normalize the data before entering it into hic_glm")
    }
    if ((is.na(contrast[1]) & is.na(coef[1])) || (!is.na(contrast[1]) & 
        !is.na(coef[1]))) {
        stop("You must enter a value for contrast or a coef, but not both")
    }
    table_list <- split(hic_table(hicexp2), hic_table(hicexp2)$chr)
    table_list <- lapply(table_list, .get_dist_tables, max.pool = max.pool)
    table_list <- do.call(c, table_list)
    dge_list <- lapply(table_list, .hictable2DGEList, covariates = meta(hicexp2))
    dge_list <- smartApply(parallel, dge_list, edgeR::estimateDisp, 
        design = design)
    fit <- smartApply(parallel, dge_list, edgeR::glmQLFit, design = design)
    if (method == "QLFTest") {
        if (is.na(coef)) {
            result <- smartApply(parallel, fit, edgeR::glmQLFTest, 
                contrast = contrast)
        }
        else {
            result <- smartApply(parallel, fit, edgeR::glmQLFTest, 
                coef = coef)
        }
    }
    if (method == "LRTest") {
        if (is.na(coef)) {
            result <- smartApply(parallel, fit, edgeR::glmLRT, 
                contrast = contrast)
        }
        else {
            result <- smartApply(parallel, fit, edgeR::glmLRT, 
                coef = coef)
        }
    }
    if (method == "Treat") {
        if (is.na(coef)) {
            result <- smartApply(parallel, fit, edgeR::glmTreat, 
                contrast = contrast, lfc = M)
        }
        else {
            result <- smartApply(parallel, fit, edgeR::glmTreat, 
                coef = coef, lfc = M)
        }
    }
    result <- mapply(.glm_reformat, result, table_list, MoreArgs = list(p.method = p.method), 
        SIMPLIFY = FALSE)
    result <- data.table::rbindlist(result)
    result <- result[order(chr, region1, region2), ]
    slot(hicexp2, "comparison") <- data.table::as.data.table(result)
    return(hicexp2)
}




data("hicexp2")
max.pool = 0.7
groups = factor(c(1,1,2,2))
hicexp2 <- fastlo(hicexp2, verbose = FALSE, parallel = FALSE)
table_list <- split(hic_table(hicexp2), hic_table(hicexp2)$chr) # split by chr

## This is the new version where we use progressive pooling of distances
# make list of tables by distance pool
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

# check to make sure number of rows for each table is not less than 500
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

table_list <- lapply(table_list, .get_dist_tables, max.pool = max.pool) # split by distance

table_list <- do.call(c, table_list)
dge_list <- lapply(table_list, .hictable2DGEList, covariates = groups)

dge_list <- lapply(dge_list, edgeR::estimateDisp, 
    design = design)
    
fit <- lapply(dge_list, edgeR::glmQLFit, design = design)

result <- lappy(fit, edgeR::glmQLFTest, 
    contrast = contrast)