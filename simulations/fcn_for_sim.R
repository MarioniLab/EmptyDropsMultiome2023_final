










SIMFUN_multiome_no_scrambling_vale <- function(raw.mat, raw.mat_atac, group1=500, group2=500, reorder.rate=0.1, down.rate.rna=0.04,  down.rate.atac=0.02) {
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
  
  
  # # Assuming all cells below k-means line are empty droplets.
  
  resampled <- raw.mat[,!above_k_means]
  
  resampled_atac <- raw.mat_atac[,!above_k_means]
  
  
  # give index names to missing colnames
  colnames(resampled)= as.character(seq(1, dim(resampled)[2]))
  colnames(resampled_atac)=as.character(seq(1, dim(resampled_atac)[2]))
  
  # Sampling real cells in group 1.
  # All things above the k-means line are assumed to be real.
  # We shuffle 10% of the genes to break any remaining "ambience".
  is.real <- which(above_k_means==1 )

  sampling1 = sample(is.real, group1, replace=TRUE)
  
  mat1 <- raw.mat[,sampling1]
  mat1 <- reorderRows(mat1, round(nrow(mat1) * reorder.rate))
  
  mat1_atac <- raw.mat_atac[,sampling1]
  mat1_atac <- reorderRows(mat1_atac, round(nrow(mat1_atac) * reorder.rate))
  
  # Sampling cells in group 2, with downsampling.
  is.good <- which(above_k_means==1 
                   & count_df_new$rna_count <14000 
                   & count_df_new$rna_count >8000
                   & count_df_new$atac_count <3000
                   & count_df_new$atac_count >800
                   )
  
  sampling2 = sample(is.good, ceiling(group2/2), replace=TRUE)
  sampling3 = sample(is.good, group2-ceiling(group2/2), replace=TRUE)
  
  mat2 <- raw.mat[,sampling2]
  mat2 <- reorderRows(mat2, round(nrow(mat2) * reorder.rate))
  mat2 <- downsampleMatrix(mat2, prop=down.rate.rna)
  
  mat2_atac <- raw.mat_atac[,sampling2]
  mat2_atac <- reorderRows(mat2_atac, round(nrow(mat2_atac) * reorder.rate))
  mat2_atac <- downsampleMatrix(mat2_atac, prop=down.rate.atac)
  
#   mat3 <- raw.mat[,sampling3]
#   mat3 <- reorderRows(mat3, round(nrow(mat3) * reorder.rate))
#   mat3 <- downsampleMatrix(mat3, prop=down.rate.rna/2)
  
#   mat3_atac <- raw.mat_atac[,sampling3]
#   mat3_atac <- reorderRows(mat3_atac, round(nrow(mat3_atac) * reorder.rate))
  
#   mat2 <- cbind(mat2, mat3)
#   mat2_atac <- cbind(mat2_atac, mat3_atac)
  
  # rename colnames to avoid overlaps
  colnames(mat1) <- paste0("g1", seq(1, dim(mat1)[2]))
  colnames(mat1_atac) <- paste0("g1", seq(1, dim(mat1_atac)[2]) )
  colnames(mat2) <- paste0("g2", seq(1, dim(mat2)[2])  )
  colnames(mat2_atac) <- paste0("g2", seq(1, dim(mat2_atac)[2]) )
  
  # Completed.
  identity=rep(0:2, c(ncol(resampled), ncol(mat1), ncol(mat2)))
  names(identity) = c(colnames(resampled), colnames(mat1), colnames(mat2))
  
  return(list(counts_rna=cbind(resampled, mat1, mat2),
              counts_atac=cbind(resampled_atac, mat1_atac, mat2_atac),
              identity=identity
  ))
  
}






SIMFUN_multiome_no_scrambling <- function(raw.mat, raw.mat_atac, group1=500, group2=500, reorder.rate=0.1, down.rate.rna=0.04,  down.rate.atac=0.02) {
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
  
  
  # # Assuming all cells below k-means line are empty droplets.
  
  resampled <- raw.mat[,!above_k_means]

  resampled_atac <- raw.mat_atac[,!above_k_means]
  
  
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
  
  return(list(counts_rna=cbind(resampled, mat1, mat2),
              counts_atac=cbind(resampled_atac, mat1_atac, mat2_atac),
              identity=identity
  ))
  
}



SIMFUN_multiome_no_scrambling_bigger_empties <- function(raw.mat, raw.mat_atac, group1=500, group2=500, reorder.rate=0.1, down.rate.rna=0.04,  down.rate.atac=0.02) {
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
  
  
  # # Assuming all cells below k-means line are empty droplets.
  
  # create ambient profiles
  ambient.prof <- rowSums(raw.mat[,!above_k_means])
  ambient.prof_atac <- rowSums(raw.mat_atac[,!above_k_means])
#   added_vector <- rmultinom(n=1, size=120, prob=ambient.prof)
#   added_vector_atac <- rmultinom(n=1, size=110, prob=ambient.prof_atac)
  
  # empties
  resampled <- raw.mat[,!above_k_means]
  resampled_atac <- raw.mat_atac[,!above_k_means]
  
#   resampled_extra = 7*resampled[, colSums(resampled)>50] 
#   resampled_atac_extra = 7*resampled_atac[, colSums(resampled_atac)>50] 
    
#   # try 2  
#   resampled_extra1 = 10*resampled[, colSums(resampled)+colSums(resampled_atac)>100 & colSums(resampled)>50] 
#   resampled_atac_extra1 = 10*resampled_atac[, colSums(resampled)+colSums(resampled_atac)>100 & colSums(resampled_atac)>50] 
#   resampled_extra2 = 7*resampled[, colSums(resampled)+colSums(resampled_atac)>100] 
#   resampled_atac_extra2 = 7*resampled_atac[, colSums(resampled)+colSums(resampled_atac)>100] 
#   resampled_extra3 = 5*resampled[, colSums(resampled)+colSums(resampled_atac)>100] 
#   resampled_atac_extra3 = 5*resampled_atac[, colSums(resampled)+colSums(resampled_atac)>100] 
#   resampled_extra4 = 3*resampled[, colSums(resampled)+colSums(resampled_atac)>100] 
#   resampled_atac_extra4 = 3*resampled_atac[, colSums(resampled)+colSums(resampled_atac)>100] 
#   resampled_extra5 = 1.5*resampled[, colSums(resampled)+colSums(resampled_atac)>100] 
#   resampled_atac_extra5 = 1.5*resampled_atac[, colSums(resampled)+colSums(resampled_atac)>100] 

#   resampled = cbind(resampled, resampled_extra1, resampled_extra2, resampled_extra3, resampled_extra4, resampled_extra5)
#   resampled_atac = cbind(resampled_atac, resampled_atac_extra1, resampled_atac_extra2, resampled_atac_extra3, resampled_atac_extra4, resampled_atac_extra5)
      
  # try 3  
#   resampled_extra = t(t(resampled[, colSums(resampled)>50 ]) * (floor(11/(150-50) *(colSums(resampled)-50)+1) ))
#   resampled_atac_extra = t(t(resampled_atac[, colSums(resampled_atac)>50])  * (floor(11/(150-50) *(colSums(resampled_atac)-50)+1) ))
#   resampled = cbind(resampled, resampled_extra)
#   resampled_atac = cbind(resampled_atac, resampled_atac_extra)
    
  # try 4 
#   resampled_extra1 = resampled[, colSums(resampled)>50 & colSums(resampled)<60]
#   resampled_extra2 = 2*resampled[, colSums(resampled)>60 & colSums(resampled)<70]
#   resampled_extra3 = 3*resampled[, colSums(resampled)>70 & colSums(resampled)<80]
#   resampled_extra4 = 4*resampled[, colSums(resampled)>80 & colSums(resampled)<90]
#   resampled_extra5 = 5*resampled[, colSums(resampled)>90 & colSums(resampled)<100]
#   resampled_extra6 = 6*resampled[, colSums(resampled)>100 & colSums(resampled)<110]
#   resampled_extra7 = 7*resampled[, colSums(resampled)>110 & colSums(resampled)<120]
#   resampled_extra8 = 8*resampled[, colSums(resampled)>120 & colSums(resampled)<130]
#   resampled_extra9 = 9*resampled[, colSums(resampled)>130 & colSums(resampled)<140]
#   resampled_extra10 = 10*resampled[, colSums(resampled)>140 ]
 
#   resampled_atac_extra1 =   resampled_atac[, colSums(resampled_atac)>50 & colSums(resampled_atac)<60]
#   resampled_atac_extra2 = 2*resampled_atac[, colSums(resampled_atac)>60 & colSums(resampled_atac)<70]
#   resampled_atac_extra3 = 3*resampled_atac[, colSums(resampled_atac)>70 & colSums(resampled_atac)<80]
#   resampled_atac_extra4 = 4*resampled_atac[, colSums(resampled_atac)>80 & colSums(resampled_atac)<90]
#   resampled_atac_extra5 = 5*resampled_atac[, colSums(resampled_atac)>90 & colSums(resampled_atac)<100]
#   resampled_atac_extra6 = 6*resampled_atac[, colSums(resampled_atac)>100 & colSums(resampled_atac)<110]
#   resampled_atac_extra7 = 7*resampled_atac[, colSums(resampled_atac)>110 & colSums(resampled_atac)<120]
#   resampled_atac_extra8 = 8*resampled_atac[, colSums(resampled_atac)>120 & colSums(resampled_atac)<130]
#   resampled_atac_extra9 = 9*resampled_atac[, colSums(resampled_atac)>130 & colSums(resampled_atac)<140]
#   resampled_atac_extra10=10*resampled_atac[, colSums(resampled_atac)>140 ]
 
#   resampled = cbind(resampled, 
#                    resampled_extra1, resampled_extra2,
#                    resampled_extra3, resampled_extra4,
#                    resampled_extra5, resampled_extra6,
#                    resampled_extra7, resampled_extra8,
#                    resampled_extra9, resampled_extra10)
#   resampled_atac = cbind(resampled_atac, 
#                    resampled_atac_extra1, resampled_atac_extra2,
#                    resampled_atac_extra3, resampled_atac_extra4,
#                    resampled_atac_extra5, resampled_atac_extra6,
#                    resampled_atac_extra7, resampled_atac_extra8,
#                    resampled_atac_extra9, resampled_atac_extra10)

  # try 5
  soup_ind = intersect(which(colSums(resampled)< 150 & colSums(resampled)>50),   which(colSums(resampled_atac)< 150 & colSums(resampled_atac)>50) )
  print("number of soup droplets: ")
  print(length(soup_ind))
  resampled_extra = resampled[ , soup_ind]
  resampled_atac_extra = resampled_atac[ , soup_ind]
  factors <- sample(2:10, length(soup_ind), replace=T)  
  factors_atac <- sample(2:40, length(soup_ind), replace=T)  
  for (i in seq(1, length(soup_ind))){
      resampled_extra[ , i]      = factors[i]*resampled_extra[ , i]
      resampled_atac_extra[ , i] = factors_atac[i]*resampled_atac_extra[ , i]
  }
  resampled = cbind(resampled, resampled_extra)
  resampled_atac = cbind(resampled_atac, resampled_atac_extra)

  
  # give index names to missing colnames
  colnames(resampled)     = as.character(seq(1, dim(resampled)[2]))
  colnames(resampled_atac)= as.character(seq(1, dim(resampled_atac)[2]))
  
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
  
  return(list(counts_rna=cbind(resampled, mat1, mat2),
              counts_atac=cbind(resampled_atac, mat1_atac, mat2_atac),
              identity=identity
  ))
  
}






SIMFUN_multiome_v3 <- function(raw.mat, raw.mat_atac, group1=500, group2=500, reorder.rate=0.1, down.rate.rna=0.04,  down.rate.atac=0.02, lower_rna=100, lower_atac=100) {
  # this version has a bimodal distribution of small cells: small in both RNA and ATAC, and very small but only in RNA
  # it scrambles the the RNA and ATAC below lower and lower_atac, but not the the rest of the droplets below the k-means line
  
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

  
  # # Assuming all cells below k-means line are empty droplets.
  
  # ambient.prof <- rowSums(raw.mat[,!above_k_means])
  # totals <- count_df_new$rna_count
  # empty.totals <- totals[!above_k_means]
  # 
  # ambient.prof_atac <- rowSums(raw.mat_atac[,!above_k_means])
  # totals_atac <- count_df_new$atac_count
  # empty.totals_atac <- totals_atac[!above_k_means]
  # 
  # # Scrambling to generate resampled empty droplets.
  # gene.ids <- rep(seq_along(ambient.prof), ambient.prof)
  # gene.ids <- sample(gene.ids)
  # cell.ids <- rep(seq_along(empty.totals), empty.totals)
  # resampled <- makeCountMatrix(gene.ids, cell.ids, all.genes=rownames(raw.mat))
  # 
  # region.ids_atac <- rep(seq_along(ambient.prof_atac), ambient.prof_atac)
  # region.ids_atac <- sample(region.ids_atac)
  # cell.ids_atac <- rep(seq_along(empty.totals_atac), empty.totals_atac)
  # resampled_atac <- makeCountMatrix(region.ids_atac, cell.ids_atac, all.genes=rownames(raw.mat_atac))
  
  

  ambient.prof <- rowSums(raw.mat[, count_df_new$atac_count<=lower_atac | count_df_new$atac_count<=lower_rna])
  totals <- count_df_new$rna_count
  empty.totals <- totals[count_df_new$atac_count<=lower_atac | count_df_new$atac_count<=lower_rna]

  ambient.prof_atac <- rowSums(raw.mat_atac[,count_df_new$atac_count<=lower_atac | count_df_new$atac_count<=lower_rna ])
  totals_atac <- count_df_new$atac_count
  empty.totals_atac <- totals_atac[count_df_new$atac_count<=lower_atac | count_df_new$atac_count<=lower_rna]

  # Scrambling to generate resampled empty droplets.
  gene.ids <- rep(seq_along(ambient.prof), ambient.prof)
  gene.ids <- sample(gene.ids)
  cell.ids <- rep(seq_along(empty.totals), empty.totals)
  resampled <- makeCountMatrix(gene.ids, cell.ids, all.genes=rownames(raw.mat))

  region.ids_atac <- rep(seq_along(ambient.prof_atac), ambient.prof_atac)
  region.ids_atac <- sample(region.ids_atac)
  cell.ids_atac <- rep(seq_along(empty.totals_atac), empty.totals_atac)
  resampled_atac <- makeCountMatrix(region.ids_atac, cell.ids_atac, all.genes=rownames(raw.mat_atac))

  # don't scramble the rest of the droplets under k_means
  rest <- raw.mat[,count_df_new$atac_count>lower_atac & count_df_new$rna_count>lower_rna & !above_k_means  ]
  rest_atac <- raw.mat_atac[,count_df_new$atac_count>lower_atac & count_df_new$rna_count>lower_rna & !above_k_means  ]
  
  # put scrambled and unscrambled together
  resampled <- cbind(resampled,rest)
  resampled_atac <- cbind(resampled_atac,rest_atac)
  
  # resampled <- raw.mat[,!above_k_means]
  # 
  # resampled_atac <- raw.mat_atac[,!above_k_means]

  
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
  
  return(list(counts_rna=cbind(resampled, mat1, mat2),
              counts_atac=cbind(resampled_atac, mat1_atac, mat2_atac),
              identity=identity
  ))
  
}




















SIMFUN_multiome_v2 <- function(raw.mat, raw.mat_atac, group1=500, group2=500, reorder.rate=0.1, down.rate.rna=0.04,  down.rate.atac=0.02, lower=100) {
  # this version has a bimodal distribution of small cells: small in both RNA and ATAC, and very small but only in RNA
  # it scrambles everything below the k-means line
  
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
  
  # Assuming all cells below k-means line are empty droplets.
  ambient.prof <- rowSums(raw.mat[,!above_k_means])
  totals <- count_df_new$rna_count
  empty.totals <- totals[!above_k_means]
  
  ambient.prof_atac <- rowSums(raw.mat_atac[,!above_k_means])
  totals_atac <- count_df_new$atac_count
  empty.totals_atac <- totals_atac[!above_k_means]
  
  # Scrambling to generate resampled empty droplets.
  gene.ids <- rep(seq_along(ambient.prof), ambient.prof)
  gene.ids <- sample(gene.ids)
  cell.ids <- rep(seq_along(empty.totals), empty.totals)
  resampled <- makeCountMatrix(gene.ids, cell.ids, all.genes=rownames(raw.mat))
  
  region.ids_atac <- rep(seq_along(ambient.prof_atac), ambient.prof_atac)
  region.ids_atac <- sample(region.ids_atac)
  cell.ids_atac <- rep(seq_along(empty.totals_atac), empty.totals_atac)
  resampled_atac <- makeCountMatrix(region.ids_atac, cell.ids_atac, all.genes=rownames(raw.mat_atac))
  
  # give index names to missing colnames
  colnames(resampled)= as.character(seq(1, dim(resampled)[2]))
  colnames(resampled_atac)=as.character(seq(1, dim(resampled_atac)[2]))
  
  # Sampling real cells in group 1.
  # All things above the k-means line are assumed to be real.
  # We shuffle 10% of the genes to break any remaining "ambience".
  is.real <- which(above_k_means==1)
  sampling1 = sample(is.real, group1, replace=TRUE)
  
  mat1 <- raw.mat[,sampling1]
  mat1 <- reorderRows(mat1, round(nrow(mat1) * reorder.rate))
  
  mat1_atac <- raw.mat_atac[,sampling1]
  mat1_atac <- reorderRows(mat1_atac, round(nrow(mat1_atac) * reorder.rate))
  
  # Sampling cells in group 2, with downsampling.
  sampling2 = sample(is.real, ceiling(group2/2), replace=TRUE)
  sampling3 = sample(is.real, group2-ceiling(group2/2), replace=TRUE)
  
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
  
  return(list(counts_rna=cbind(resampled, mat1, mat2),
              counts_atac=cbind(resampled_atac, mat1_atac, mat2_atac),
              identity=identity
  ))
  
}


SIMFUN_multiome <- function(raw.mat, raw.mat_atac, group1=500, group2=500, reorder.rate=0.1, down.rate=0.1, lower=100) {
  
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

  # Assuming all cells below k-means line are empty droplets.
  ambient.prof <- rowSums(raw.mat[,!above_k_means])
  totals <- count_df_new$rna_count
  empty.totals <- totals[!above_k_means]
  
  ambient.prof_atac <- rowSums(raw.mat_atac[,!above_k_means])
  totals_atac <- count_df_new$atac_count
  empty.totals_atac <- totals_atac[!above_k_means]
  
  # Scrambling to generate resampled empty droplets.
  gene.ids <- rep(seq_along(ambient.prof), ambient.prof)
  gene.ids <- sample(gene.ids)
  cell.ids <- rep(seq_along(empty.totals), empty.totals)
  resampled <- makeCountMatrix(gene.ids, cell.ids, all.genes=rownames(raw.mat))
  
  region.ids_atac <- rep(seq_along(ambient.prof_atac), ambient.prof_atac)
  region.ids_atac <- sample(region.ids_atac)
  cell.ids_atac <- rep(seq_along(empty.totals_atac), empty.totals_atac)
  resampled_atac <- makeCountMatrix(region.ids_atac, cell.ids_atac, all.genes=rownames(raw.mat_atac))
  
  # give index names to missing colnames
  colnames(resampled)= as.character(seq(1, dim(resampled)[2]))
  colnames(resampled_atac)=as.character(seq(1, dim(resampled_atac)[2]))
  
  # Sampling real cells in group 1.
  # All things above the k-means line are assumed to be real.
  # We shuffle 10% of the genes to break any remaining "ambience".
  is.real <- which(above_k_means==1)
  sampling1 = sample(is.real, group1, replace=TRUE)
  
  mat1 <- raw.mat[,sampling1]
  mat1 <- reorderRows(mat1, round(nrow(mat1) * reorder.rate))
  
  mat1_atac <- raw.mat_atac[,sampling1]
  mat1_atac <- reorderRows(mat1_atac, round(nrow(mat1_atac) * reorder.rate))
  
  # Sampling cells in group 2, with downsampling.
  sampling2 = sample(is.real, group2, replace=TRUE)
  
  mat2 <- raw.mat[,sampling2]
  mat2 <- reorderRows(mat2, round(nrow(mat2) * reorder.rate))
  mat2 <- downsampleMatrix(mat2, prop=down.rate*2)
  
  mat2_atac <- raw.mat_atac[,sampling2]
  mat2_atac <- reorderRows(mat2_atac, round(nrow(mat2_atac) * reorder.rate))
  mat2_atac <- downsampleMatrix(mat2_atac, prop=down.rate)
  
  # rename colnames to avoid overlaps
  colnames(mat1) <- paste0("g1", seq(1, dim(mat1)[2]))
  colnames(mat1_atac) <- paste0("g1", seq(1, dim(mat1_atac)[2]) )
  colnames(mat2) <- paste0("g2", seq(1, dim(mat2)[2])  )
  colnames(mat2_atac) <- paste0("g2", seq(1, dim(mat2_atac)[2]) )
  
  # Completed.
  identity=rep(0:2, c(ncol(resampled), ncol(mat1), ncol(mat2)))
  names(identity) = c(colnames(resampled), colnames(mat1), colnames(mat2))

  return(list(counts_rna=cbind(resampled, mat1, mat2),
              counts_atac=cbind(resampled_atac, mat1_atac, mat2_atac),
              identity=identity
              ))
  
}




SIMFUN <- function(raw.mat, group1=500, group2=500, reorder.rate=0.1, down.rate=0.1, lower=100) {
  stats <- barcodeRanks(raw.mat)
  
  # Assuming all cells below the inflection point are empty droplets.
  totals <- stats$total
  empty.limit <- stats@metadata[["inflection"]]
  ambient.prof <- rowSums(raw.mat[,totals <= empty.limit])
  empty.totals <- totals[totals <= empty.limit]
  
  # Scrambling to generate resampled empty droplets.
  gene.ids <- rep(seq_along(ambient.prof), ambient.prof)
  gene.ids <- sample(gene.ids)
  cell.ids <- rep(seq_along(empty.totals), empty.totals)
  resampled <- makeCountMatrix(gene.ids, cell.ids, all.genes=rownames(raw.mat))
  
  # Sampling real cells in group 1.
  # All things above the inflection point are assumed to be real.
  # We shuffle 10% of the genes to break any remaining "ambience".
  is.real <- which(empty.limit <= totals)
  #print("debug1")
  mat1 <- raw.mat[,sample(is.real, group1, replace=TRUE)]
  #print("debug2")
  mat1 <- reorderRows(mat1, round(nrow(mat1) * reorder.rate))
  
  # Sampling cells in group 2, with downsampling.
  mat2 <- raw.mat[,sample(is.real, group2, replace=TRUE)]
  mat2 <- reorderRows(mat2, round(nrow(mat2) * reorder.rate))
  mat2 <- downsampleMatrix(mat2, prop=down.rate)
  
  # Completed.
  return(list(counts=cbind(resampled, mat1, mat2),
              identity=rep(0:2, c(ncol(resampled), ncol(mat1), ncol(mat2)))))
}




















# 
# 
# raw.mat <- sce[["Gene Expression"]]
# raw.mat_atac <- sce[["Peaks"]]
# 
# 
# stats <- barcodeRanks(raw.mat)
# stats_atac <- barcodeRanks(raw.mat_atac)
# # 
# p1 <- plot(stats$rank, stats$total, log="xy", xlab="Rank", ylab="Total")
# p2 <- plot(stats_atac$rank, stats_atac$total, log="xy", xlab="Rank", ylab="Total")
# 
# p1 | p2
# 
# # Assuming all cells below the inflection point are empty droplets.
# totals <- stats$total
# empty.limit <- stats@metadata[["inflection"]]
# ambient.prof <- rowSums(raw.mat[,totals <= empty.limit])
# empty.totals <- totals[totals <= empty.limit]
# 
# # Scrambling to generate resampled empty droplets.
# gene.ids <- rep(seq_along(ambient.prof), ambient.prof)
# gene.ids <- sample(gene.ids)  # randomly permute gene.ids
# cell.ids <- rep(seq_along(empty.totals), empty.totals)
# resampled <- makeCountMatrix(gene.ids, cell.ids, all.genes=rownames(raw.mat))
# 
# # Sampling real cells in group 1.
# # All things above the inflection point are assumed to be real.
# # We shuffle 10% of the genes to break any remaining "ambience".
# is.real <- which(empty.limit <= totals)
# mat1 <- raw.mat[,sample(is.real, group1, replace=TRUE)]
# mat1 <- reorderRows(mat1, round(nrow(mat1) * reorder.rate))
# 
# # Sampling cells in group 2, with downsampling.
# mat2 <- raw.mat[,sample(is.real, group2, replace=TRUE)]
# mat2 <- reorderRows(mat2, round(nrow(mat2) * reorder.rate))
# mat2 <- downsampleMatrix(mat2, prop=down.rate)
# 
# # Completed.
# return(list(counts=cbind(resampled, mat1, mat2),
#             identity=rep(0:2, c(ncol(resampled), ncol(mat1), ncol(mat2)))))
# 

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
