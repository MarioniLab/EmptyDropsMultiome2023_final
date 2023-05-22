setwd("/mnt/beegfs6/home3/ahringer/em613/analysis/multiomics/emptryDrops_multiome2023")
library(eDv2)
library(DropletUtils)
library(Matrix)
library(rstudioapi)
library(Seurat)
library(mixtools)
library(tidyverse)
source("simulations/fcn_for_sim.R")
# library("optparse")
library(Seurat)
library(Signac)
library(stringr)


sce <- Read10X_h5(file.path("data/input", "pbmc_gran_sorted/raw_feature_bc_matrix.h5"))
count_matrix_atac <- sce[["Peaks"]]


# calculate sizes of peaks and order by them
ends_all <- str_split(  lapply(str_split(rownames(  count_matrix_atac  ), pattern = ":"), `[[`, 2)  ,  pattern = "-" )
peaks_sizes = unlist(lapply(ends_all, FUN = function(x) strtoi(x[2])-strtoi(x[1])  )   )
op = order(peaks_sizes)
count_matrix_atac = count_matrix_atac[op,]
n_shortest = 14000

# calculate library sizes of barcodes and order by them
totals = colSums(count_matrix_atac)
ot = order(totals)
count_matrix_atac = count_matrix_atac[,ot]

# create simulated matrix
real_cells_cols = sample( which(  colSums(count_matrix_atac)>20000 & colSums(count_matrix_atac)<35000  ), 2000  )
sim_matrix = count_matrix_atac[,real_cells_cols]
sim_matrix[seq(1,n_shortest) , seq(1001,2000) ] = 0
sim_matrix[seq(n_shortest+1,2*n_shortest) , seq(1,1000)] = 0

empties_cols = sample( which(  colSums(count_matrix_atac)<200  )  )
empt_matrix = count_matrix_atac[,empties_cols]
half = round( length(empties_cols)/2 )
empt_matrix[seq(1,n_shortest), seq(half, length(empties_cols))] = 0
empt_matrix[seq(n_shortest+1,2*n_shortest), seq(half, length(empties_cols))] = 0

# create nonpeak column
meta <- read.table(file.path("data/input", "pbmc_gran_sorted/pbmc_granulocyte_sorted_10k_per_barcode_metrics.csv"), header=T, sep=",")
nonpeak = (2*meta$atac_fragments - meta$atac_peak_region_cutsites)/(2*meta$atac_fragments)

max_frip_of_excluded = max(nonpeak[meta$excluded_reason==2]  )
min_frip_of_cells = min(nonpeak[meta$is_cell==1]  )
hist(nonpeak, ylim= c(0,10000), breaks=500)

CSiP = (meta$atac_peak_region_cutsites)/(2*meta$atac_fragments)
max_CSiP_of_excluded = max(CSiP[meta$excluded_reason==2]  )
max_CSiP_of_excluded
min_CSiP_of_cells = min(CSiP[meta$is_cell==1]  )
min_CSiP_of_cells
hist(nonpeak, ylim= c(0,10000), breaks=500)


nonpeak = meta$atac_fragments - meta$atac_peak_region_fragments
names(nonpeak) = meta$barcode
nonpeak = nonpeak[colnames(count_matrix_atac)]
count_matrix_w_nonpeak = rbind(count_matrix_atac, nonpeak)



