SIMFUN_multiome_no_scrambling_more_empties  <- function(raw.mat, raw.mat_atac, group1=500, group2=500, reorder.rate=0.1, down.rate.rna=0.04,  down.rate.atac=0.02) {
  # this version has a bimodal distribution of small cells: small in both RNA and ATAC, and very small but only in RNA
  
  # remove NA values from the count matrices
  raw.mat[ is.na(raw.mat)  ] = 0
  raw.mat_atac[ is.na(raw.mat_atac)  ] = 0
  
  # change column names to the union of column names
  outlist = match_colnames(raw.mat, raw.mat_atac)
  raw.mat = outlist$mat1_new
  raw.mat_atac = outlist$mat2_new
  
  # order column names alphabetically in both count matrices
  raw.mat = raw.mat[,order(colnames(raw.mat))]
  raw.mat_atac = raw.mat_atac[,colnames(raw.mat)]
  
  # create dataframe to call_cells by k-means
  count_df <- data.frame("rownames" = colnames(raw.mat),
                         "atac_count" = unname(colSums(raw.mat_atac)),
                         "rna_count" = unname(colSums(raw.mat)),
                         "excluded"  = rep(F, length(colSums(raw.mat))),
                         "is_cell"  =  rep(0, length(colSums(raw.mat)))
  )
  list_output <- call_cells(count_df)
  count_df_new <- list_output[[1]]
  count_df_new <- count_df_new[order(count_df_new$rownames),]
  above_k_means <- count_df_new$is_cell
  k_means_eq <- list_output[[2]]
  print(k_means_eq)  
  
  
  
  
  # # Assuming all cells below k-means line are empty droplets.
  
  # create more empties
  ambient.prof <- 4*rowSums(raw.mat[,colSums(raw.mat)<=150 & colSums(raw.mat)>=10])
  totals <- count_df_new$rna_count
  
  ambient.prof_atac <- 4*rowSums(raw.mat_atac[,colSums(raw.mat_atac)<=90])
  totals_atac <- count_df_new$atac_count
  
  print("simulate empties")
  max_rna = log10(max(totals)/2)
  min_rna =  log10(min(totals)+20)
  max_atac = log10(max(totals_atac)/10)
  min_atac =  log10(min(totals_atac)+20)
  # new_rna_cts = floor(runif(400000, min=min_rna, max=max_rna))
  # new_atac_cts = floor(runif(400000, min=min_atac, max=max_atac))
  new_rna_cts = floor(10^(runif(400000, min=min_rna, max=max_rna)) )
  new_atac_cts = floor(10^(runif(400000, min=min_atac, max=max_atac)))
  
  keep = which(log10(new_rna_cts+0.1) - k_means_eq[2] * log10(new_atac_cts+0.1 )  > k_means_eq[1]-1.5 & 
                 log10(new_rna_cts+0.1) - k_means_eq[2] * log10(new_atac_cts+0.1 )  < k_means_eq[1] + 3 &
                 log10(new_rna_cts+0.1) + 1/k_means_eq[2] * log10(new_atac_cts+0.1 )  <  1.1 ) 
  
  keep = keep[seq(1,4000)]
  new_rna_cts = new_rna_cts[keep]
  new_atac_cts = new_atac_cts[keep]
  
  
  # slow way of rescaling amb prof
  # regMat <- matrix(, nrow = length(ambient.prof), ncol = 0)
  # regMat_atac <- matrix(, nrow = length(ambient.prof_atac), ncol = 0)
  # 
  # start_time <- Sys.time()
  # 
  # for (i in seq_along(new_rna_cts)){
  #   #print(i)
  #   new_empt <- floor( ambient.prof * new_rna_cts[i] / sum(ambient.prof) )
  #   regMat <- cbind(regMat, new_empt)
  #   
  #   new_empt_atac <- floor( ambient.prof_atac * new_atac_cts[i] / sum(ambient.prof_atac) )
  #   regMat_atac <- cbind(regMat_atac, new_empt_atac)
  # }
  # 
  # end_time <- Sys.time()
  # print(end_time - start_time)
  # start_time <- Sys.time()
  # 
  # resampled <- Matrix(regMat, sparse = TRUE)
  # resampled_atac <- Matrix(regMat_atac, sparse = TRUE)
  # 
  
  # fast way of rescaling amb prof
  start_time <- Sys.time()
  listMat <- lapply(new_rna_cts, function(x) ambient.prof * x / sum(ambient.prof) )
  regMat <- do.call(cbind, listMat)
  listMat_atac <- lapply(new_atac_cts, function(x) ambient.prof_atac * x / sum(ambient.prof_atac) )
  regMat_atac <- do.call(cbind, listMat_atac)
  resampled <- Matrix(regMat, sparse = TRUE)
  resampled_atac <- Matrix(regMat_atac, sparse = TRUE)
  end_time <- Sys.time()
  print(end_time - start_time)
  
  
  # sample with multinomial distribution
  # gene.ids <- rep(seq_along(ambient.prof), ambient.prof)
  # gene.ids <- sample(gene.ids)
  # cell.ids <- rep(seq_along(new_rna_cts), new_rna_cts)
  # gene.ids <- gene.ids[seq_along(cell.ids)]
  # resampled <- makeCountMatrix(gene.ids, cell.ids, all.genes=rownames(raw.mat))
  # 
  # gene.ids <- rep(seq_along(ambient.prof_atac), ambient.prof_atac)
  # gene.ids <- sample(gene.ids)
  # cell.ids <- rep(seq_along(new_atac_cts), new_atac_cts)
  # gene.ids <- gene.ids[seq_along(cell.ids)]
  # resampled_atac <- makeCountMatrix(gene.ids, cell.ids, all.genes=rownames(raw.mat_atac))
  # 
  
  
  plot(log10(new_atac_cts + 0.1), log10(new_rna_cts + 0.1),
       #     xlim=c(0,1510) , ylim=c(0,8),
       pch=".",
       # cex=2,
       xlab="log10(atac_count)", ylab="log10(rna_count)",
       main=paste0("simulated empties")
  )+theme(axis.title.x=element_text(size=14,face="bold"),
          axis.title.y=element_text(size=14,face="bold"))
  
  
  resampled <- cbind(raw.mat[,!above_k_means], resampled)
  resampled_atac <- cbind( raw.mat_atac[,!above_k_means], resampled_atac)
  
  print("finished simulating empties")
  
  # give index names to missing colnames
  colnames(resampled)= as.character(seq(1, dim(resampled)[2]))
  colnames(resampled_atac)=as.character(seq(1, dim(resampled_atac)[2]))
  
  # Sampling real cells in group 1.
  # All things above the k-means line are assumed to be real.
  # We shuffle 10% of the genes to break any remaining "ambience".
  is.real <- which(above_k_means==1)
  sampling1 = sample(is.real, group1, replace=FALSE)
  
  mat1 <- raw.mat[,sampling1]
  mat1 <- reorderRows(mat1, round(nrow(mat1) * reorder.rate))
  
  mat1_atac <- raw.mat_atac[,sampling1]
  mat1_atac <- reorderRows(mat1_atac, round(nrow(mat1_atac) * reorder.rate))
  
  # Sampling cells in group 2, with downsampling.
  sampling2 = sample(is.real, ceiling(group2/2), replace=FALSE)
  sampling3 = sample(is.real, group2-ceiling(group2/2), replace=FALSE)
  
  mat2 <- raw.mat[,sampling2]
  mat2 <- reorderRows(mat2, round(nrow(mat2) * reorder.rate))
  mat2 <- downsampleMatrix(mat2, prop=down.rate.rna)
  
  mat2_atac <- raw.mat_atac[,sampling2]
  mat2_atac <- reorderRows(mat2_atac, round(nrow(mat2_atac) * reorder.rate))
  mat2_atac <- downsampleMatrix(mat2_atac, prop=down.rate.atac)
  
  mat3 <- raw.mat[,sampling3]
  mat3 <- reorderRows(mat3, round(nrow(mat3) * reorder.rate))
  mat3 <- downsampleMatrix(mat3, prop=down.rate.rna/2)
  
  mat3_atac <- raw.mat_atac[,sampling3]
  mat3_atac <- reorderRows(mat3_atac, round(nrow(mat3_atac) * reorder.rate))
  
  mat2 <- cbind(mat2, mat3)
  mat2_atac <- cbind(mat2_atac, mat3_atac)
  
  # rename colnames to avoid overlaps
  colnames(mat1) <- paste0("g1", seq(1, dim(mat1)[2]))
  colnames(mat1_atac) <- paste0("g1", seq(1, dim(mat1_atac)[2]) )
  colnames(mat2) <- paste0("g2", seq(1, dim(mat2)[2])  )
  colnames(mat2_atac) <- paste0("g2", seq(1, dim(mat2_atac)[2]) )
  
  # Completed.
  identity=rep(0:2, c(ncol(resampled), ncol(mat1), ncol(mat2)))
  names(identity) = c(colnames(resampled), colnames(mat1), colnames(mat2))
  
  end_time <- Sys.time()
  print(end_time - start_time)
  
  plot(log10(colSums(cbind(resampled_atac, mat1_atac, mat2_atac)) + 0.1), log10(colSums(cbind(resampled, mat1, mat2)) + 0.1),
       #     xlim=c(0,1510) , ylim=c(0,8),
       pch=".",
       # cex=2,
       xlab="log10(atac_count)", ylab="log10(rna_count)",
       main=paste0("simulated droplets")
  )+theme(axis.title.x=element_text(size=14,face="bold"),
          axis.title.y=element_text(size=14,face="bold"))
  
  
  
  return(list(counts_rna=cbind(resampled, mat1, mat2),
              counts_atac=cbind(resampled_atac, mat1_atac, mat2_atac),
              identity=identity
  ))
  
}













SIMFUN_multiome_no_scrambling_vale_more_empties  <- function(raw.mat, raw.mat_atac, group1=500, group2=500, reorder.rate=0.1, down.rate.rna=0.04,  down.rate.atac=0.02) {
  # this version has a bimodal distribution of small cells: small in both RNA and ATAC, and very small but only in RNA
  
  # remove NA values from the count matrices
  raw.mat[ is.na(raw.mat)  ] = 0
  raw.mat_atac[ is.na(raw.mat_atac)  ] = 0
  
  # change column names to the union of column names
  outlist = match_colnames(raw.mat, raw.mat_atac)
  raw.mat = outlist$mat1_new
  raw.mat_atac = outlist$mat2_new
  
  # order column names alphabetically in both count matrices
  raw.mat = raw.mat[,order(colnames(raw.mat))]
  raw.mat_atac = raw.mat_atac[,colnames(raw.mat)]
  
  # create dataframe to call_cells by k-means
  count_df <- data.frame("rownames" = colnames(raw.mat),
                         "atac_count" = unname(colSums(raw.mat_atac)),
                         "rna_count" = unname(colSums(raw.mat)),
                         "excluded"  = rep(F, length(colSums(raw.mat))),
                         "is_cell"  =  rep(0, length(colSums(raw.mat)))
  )
  list_output <- call_cells(count_df)
  count_df_new <- list_output[[1]]
  count_df_new <- count_df_new[order(count_df_new$rownames),]
  above_k_means <- count_df_new$is_cell
  k_means_eq <- list_output[[2]]
  print(k_means_eq)  
  
  
  
  
  # # Assuming all cells below k-means line are empty droplets.
  
  # create more empties
  ambient.prof <- 4*rowSums(raw.mat[,!above_k_means & colSums(raw.mat)<1000])
  totals <- count_df_new$rna_count
  empty.totals <- totals[!above_k_means]
  
  ambient.prof_atac <- 4*rowSums(raw.mat_atac[,!above_k_means])
  totals_atac <- count_df_new$atac_count
  empty.totals_atac <- totals_atac[!above_k_means]
  
  
  
  
  print("simulate empties")
  max_rna = log10(max(totals)/20)
  min_rna =  log10(min(totals)+20)
  max_atac = log10(max(totals_atac)/20)
  min_atac =  log10(min(totals_atac)+20)
  # new_rna_cts = floor(runif(400000, min=min_rna, max=max_rna))
  # new_atac_cts = floor(runif(400000, min=min_atac, max=max_atac))
  new_rna_cts = floor(10^(runif(400000, min=min_rna, max=max_rna)) )
  new_atac_cts = floor(10^(runif(400000, min=min_atac, max=max_atac)))
  
  keep = which(log10(new_rna_cts+0.1) - k_means_eq[2] * log10(new_atac_cts+0.1 )  > k_means_eq[1]-1.5 & 
                 log10(new_rna_cts+0.1) - k_means_eq[2] * log10(new_atac_cts+0.1 )  < k_means_eq[1] + 2 )#&
  #log10(new_rna_cts+0.1) + 1/k_means_eq[2] * log10(new_atac_cts+0.1 )  <  1 ) 
  
  keep = keep[seq(1,4000)]
  new_rna_cts = new_rna_cts[keep]
  new_atac_cts = new_atac_cts[keep]
  
  
  # fast way of rescaling amb prof
  # start_time <- Sys.time()
  # listMat <- lapply(new_rna_cts, function(x) ambient.prof * x / sum(ambient.prof) )
  # regMat <- do.call(cbind, listMat)
  # listMat_atac <- lapply(new_atac_cts, function(x) ambient.prof_atac * x / sum(ambient.prof_atac) )
  # regMat_atac <- do.call(cbind, listMat_atac)
  # resampled <- Matrix(regMat, sparse = TRUE)
  # resampled_atac <- Matrix(regMat_atac, sparse = TRUE)
  # end_time <- Sys.time()
  # print(end_time - start_time)  
  
  
  # generate empties by multinomial sampling  
  start_time <- Sys.time()
  gene.ids <- rep(seq_along(ambient.prof), ambient.prof)
  gene.ids <- sample(gene.ids)
  cell.ids <- rep(seq_along(new_rna_cts), new_rna_cts)
  gene.ids <- gene.ids[seq_along(cell.ids)]
  resampled <- makeCountMatrix(gene.ids, cell.ids, all.genes=rownames(raw.mat))
  
  gene.ids <- rep(seq_along(ambient.prof_atac), ambient.prof_atac)
  gene.ids <- sample(gene.ids)
  cell.ids <- rep(seq_along(new_atac_cts), new_atac_cts)
  gene.ids <- gene.ids[seq_along(cell.ids)]
  resampled_atac <- makeCountMatrix(gene.ids, cell.ids, all.genes=rownames(raw.mat_atac))
  end_time <- Sys.time()
  print(end_time - start_time)
  
  
  plot(log10(new_atac_cts + 0.1), log10(new_rna_cts + 0.1),
       #     xlim=c(0,1510) , ylim=c(0,8),
       pch=".",
       # cex=2,
       xlab="log10(atac_count)", ylab="log10(rna_count)",
       main=paste0("simulated empties")
  )+theme(axis.title.x=element_text(size=14,face="bold"),
          axis.title.y=element_text(size=14,face="bold"))
  
  
  resampled <- cbind(raw.mat[,!above_k_means], resampled)
  resampled_atac <- cbind( raw.mat_atac[,!above_k_means], resampled_atac)
  
  print("finished simulating empties")
  
  # give index names to missing colnames
  colnames(resampled)= as.character(seq(1, dim(resampled)[2]))
  colnames(resampled_atac)=as.character(seq(1, dim(resampled_atac)[2]))
  
  # Sampling real cells in group 1.
  # All things above the k-means line are assumed to be real.
  # We shuffle 10% of the genes to break any remaining "ambience".
  is.real <- which(above_k_means==1)
  sampling1 = sample(is.real, group1, replace=FALSE)
  
  mat1 <- raw.mat[,sampling1]
  mat1 <- reorderRows(mat1, round(nrow(mat1) * reorder.rate))
  
  mat1_atac <- raw.mat_atac[,sampling1]
  mat1_atac <- reorderRows(mat1_atac, round(nrow(mat1_atac) * reorder.rate))
  
  # Sampling cells in group 2, with downsampling.
#   is.good <- which(above_k_means==1 
#                    & count_df_new$rna_count <14000 
#                    & count_df_new$rna_count >8000
#                    & count_df_new$atac_count <3000
#                    & count_df_new$atac_count >800
#   )
  
  sampling2 = sample(is.real, ceiling(group2/2), replace=TRUE)
  sampling3 = sample(is.real, group2-ceiling(group2/2), replace=TRUE)
  
  # sampling2 = sample(is.real, ceiling(group2/2), replace=FALSE)
  # sampling3 = sample(is.real, group2-ceiling(group2/2), replace=FALSE)
  # 
  mat2 <- raw.mat[,sampling2]
  mat2 <- reorderRows(mat2, round(nrow(mat2) * reorder.rate))
  mat2 <- downsampleMatrix(mat2, prop=down.rate.rna)
  
  mat2_atac <- raw.mat_atac[,sampling2]
  mat2_atac <- reorderRows(mat2_atac, round(nrow(mat2_atac) * reorder.rate))
  mat2_atac <- downsampleMatrix(mat2_atac, prop=down.rate.atac)
  
  mat3 <- raw.mat[,sampling3]
  mat3 <- reorderRows(mat3, round(nrow(mat3) * reorder.rate))
  mat3 <- downsampleMatrix(mat3, prop=down.rate.rna/2)
  
  mat3_atac <- raw.mat_atac[,sampling3]
  mat3_atac <- reorderRows(mat3_atac, round(nrow(mat3_atac) * reorder.rate))
  
  mat2 <- cbind(mat2, mat3)
  mat2_atac <- cbind(mat2_atac, mat3_atac)
  
  # rename colnames to avoid overlaps
  colnames(mat1) <- paste0("g1", seq(1, dim(mat1)[2]))
  colnames(mat1_atac) <- paste0("g1", seq(1, dim(mat1_atac)[2]) )
  colnames(mat2) <- paste0("g2", seq(1, dim(mat2)[2])  )
  colnames(mat2_atac) <- paste0("g2", seq(1, dim(mat2_atac)[2]) )
  
  
  plot(log10( colSums(cbind(resampled_atac, mat1_atac, mat2_atac)) + 0.1), log10(  colSums(cbind(resampled, mat1, mat2)) + 0.1),
       #     xlim=c(0,1510) , ylim=c(0,8),
       pch=".",
       # cex=2,
       col=factor(rep(0:2, c(ncol(resampled), ncol(mat1), ncol(mat2))) ),
       xlab="log10(atac_count)", ylab="log10(rna_count)",
       main=paste0("simulated empties")
  )+theme(axis.title.x=element_text(size=14,face="bold"),
          axis.title.y=element_text(size=14,face="bold"))
  
  
  # Completed.
  identity=rep(0:2, c(ncol(resampled), ncol(mat1), ncol(mat2)))
  names(identity) = c(colnames(resampled), colnames(mat1), colnames(mat2))
  
  return(list(counts_rna=cbind(resampled, mat1, mat2),
              counts_atac=cbind(resampled_atac, mat1_atac, mat2_atac),
              identity=identity
  ))
  
}




















match_colnames <- function(mat1, mat2){
  
  colnames1 <- colnames(mat1)
  colnames2 <- colnames(mat2)
  union <- union(colnames1,colnames2)
  colnames1_minus_2 <- setdiff(colnames1, colnames2)
  colnames2_minus_1 <- setdiff(colnames2, colnames1)
  
  more_mat1 <- matrix(0, dim(mat1)[1] ,length(colnames2_minus_1))
  mat1_new <- cbind(mat1, more_mat1)
  
  more_mat2 <- matrix(0, dim(mat2)[1] ,length(colnames1_minus_2))
  mat2_new <- cbind(mat2, more_mat2)
  
  return(list("mat1_new"=mat1_new, "mat2_new"=mat2_new )  )
}


reorderRows <- function(mat, N) {
  i <- seq_len(nrow(mat))
  choice <- sample(nrow(mat), N)
  i[choice] <- sample(choice)
  mat[i,,drop=FALSE]
}

plotBarcodes <- function(ranks, totals, fitted=NULL, subset=NULL, ...) {
  xlim <- range(ranks[ranks>0])
  ylim <- range(totals[totals>0])
  
  # Dropping non-unique points to save space.
  # Subsetting performed after range() to create comparable plots, if desired.
  keep <- !duplicated(totals)
  if (!is.null(subset)) {
    alt.keep <- keep
    keep[] <- FALSE
    keep[subset] <- alt.keep[subset] 
  }
  Rx <- ranks[keep]
  Tx <- totals[keep]
  
  # Actually making the plot, also plotting the fitted values if requested.
  plot(Rx, Tx, log="xy", xlab="Rank", ylab="Total count", xlim=xlim, ylim=ylim, ...)
  if (!is.null(fitted)) {  
    Fx <- fitted[keep]
    o <- order(Rx)
    lines(Rx[o], Fx[o], col="red", lwd=2)
  }
  return(invisible(NULL))
}

assessMethod <- function(keep, identity) {
  g1 <- sum(keep & identity==1, na.rm=TRUE)/sum(identity==1)
  g2 <- sum(keep & identity==2, na.rm=TRUE)/sum(identity==2)
  fdr <- sum(identity==0 & keep, na.rm=TRUE)/sum(keep, na.rm=TRUE)
  return(c(G1=g1, G2=g2, FDR=fdr))
}

createRocPts <- function(stats, identity, n=20000) {
  o <- order(stats)
  stats <- stats[o]
  identity <- factor(identity[o])
  
  chunked <- cut(stats, n)
  by.chunk <- split(identity, chunked)
  N <- do.call(rbind, lapply(by.chunk, table))
  
  apply(N, 2, cumsum)
}

plotHistogramOutline <- function(breaks, heights, ...) {
  x <- rep(breaks, each=2)
  y <- c(0, rep(heights, each=2), 0)
  lines(x, y, ...)
}

#colors <- c("EmptyDrops"="salmon", "CellRanger"="dodgerblue", "Knee point"="orange")
colors <- c("eD-multiome"="salmon", "eD-rna"="green", "eD-atac"="orange", "CellRanger-arc"="dodgerblue")
