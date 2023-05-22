setwd("/mnt/beegfs6/home3/ahringer/em613/analysis/multiomics/emptryDrops_multiome2023_eDv3")
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

current_date= paste(unlist(strsplit(as.character(Sys.Date()), "-")), collapse="")
opath <- paste0("data/output/sim_frip/", current_date, "/results-sim")
old_date = "20230310"

dir.create(opath,recursive=TRUE)

set.seed(73953024)

ALLFILES <- c(#"pbmc_gran_sorted/raw_feature_bc_matrix.h5",
  # "gastr_d3.5_A/gastr_d3_5_multiome_A_L001_raw_feature_bc_matrix.h5",
  # "gastr_d3.5_B/gastr_d3_5_multiome_B_L001_raw_feature_bc_matrix.h5",
  # "gastr_d4.5_A/gastr_d4_5_multiome_A_L001_raw_feature_bc_matrix.h5",
  # "gastr_d4.5_B/gastr_d4_5_multiome_B_L001_raw_feature_bc_matrix.h5",
  # "gastr_d5_A/gastr_d5_multiome_A_L001_raw_feature_bc_matrix.h5",
  # "gastr_d5_B/gastr_d5_multiome_B_L001_raw_feature_bc_matrix.h5",
  # "gastr_d4_A/raw_feature_bc_matrix.h5",
  # "gastr_d4_B/gastr_d4_multiome_B_L001_raw_feature_bc_matrix.h5"
  "valentina_8176/FCA_GND10288176_raw_feature_bc_matrix.h5",
  "valentina_8177/FCA_GND10288177_raw_feature_bc_matrix.h5",
  "valentina_8178/FCA_GND10288178_raw_feature_bc_matrix.h5",
  "valentina_8179/FCA_GND10288179_raw_feature_bc_matrix.h5",
  "valentina_8180/FCA_GND10288180_raw_feature_bc_matrix.h5"

)



METADATA <- c(
  "valentina_8176/FCA_GND10288176_per_barcode_metrics.csv",
  "valentina_8177/FCA_GND10288177_per_barcode_metrics.csv",
  "valentina_8178/FCA_GND10288178_per_barcode_metrics.csv",
  "valentina_8179/FCA_GND10288179_per_barcode_metrics.csv",
  "valentina_8180/FCA_GND10288180_per_barcode_metrics.csv"
  
)






for (i in seq_along(ALLFILES) ) { 
  
  fname = ALLFILES[i]
  fmeta = METADATA[i]
  
  sce <- Read10X_h5( file.path("data/input", fname) )
  meta <- read.table( file.path("data/input", fmeta) , header = T, sep=",")
  
  stub <- sub("/.*", "", fname, "_qc")
  
  ffile <- file.path(opath, paste0(stub, "_frip.pdf"))  

  count_matrix_atac <- sce[["Peaks"]]
  
  nonpeak = meta$atac_fragments - meta$atac_peak_region_fragments
  names(nonpeak) = meta$barcode
  nonpeak = nonpeak[colnames(count_matrix_atac)]
  
  cutsites_in_peaks = meta$atac_peak_region_cutsites
  names(cutsites_in_peaks) = meta$barcode
  cutsites_in_peaks = cutsites_in_peaks[colnames(count_matrix_atac)]
  
  count_matrix_w_nonpeak = rbind(count_matrix_atac, nonpeak)
  
  small_peaks = sort(rowSums(count_matrix_atac))[seq(1, round(dim(count_matrix_atac)[1]*0.1)  )]
  
  # calculate sizes of peaks
  ends = str_split(  lapply(str_split(names(small_peaks), pattern = ":"), `[[`, 2)  ,  pattern = "-" )
  ends_all <- str_split(  lapply(str_split(rownames(  count_matrix_atac  ), pattern = ":"), `[[`, 2)  ,  pattern = "-" )
  # calculate fraction of genome in peaks with and without padding
  GenInPeaks_no_pad <- sum(  unlist(lapply(ends_all, FUN = function(x) strtoi(x[2])-strtoi(x[1])  ) )   )*100 /(3*10^9)
  GenInPeaks_500totalpad <- sum(  unlist(lapply(ends_all, FUN = function(x) strtoi(x[2])-strtoi(x[1])+500  ) )   )*100 /(3*10^9)
  GenInPadPeaks_2000totalpad <- sum(  unlist(lapply(ends_all, FUN = function(x) strtoi(x[2])-strtoi(x[1])+2*4000  ) )   )*100 /(3*10^9)
  
  pdf(ffile)
  
  plot.new()
  text(x=.1, y=.1, paste0("GenInPeaks_no_pad=", GenInPeaks_no_pad))
  text(x=.25, y=0.6, paste0("GenInPeaks_500totalpad=", GenInPeaks_500totalpad))
  text(x=0.5, y=0.8, paste0("GenInPadPeaks_2000totalpad=", GenInPadPeaks_2000totalpad))
  
  small_peaks_sizes = unlist(lapply(ends, FUN = function(x) strtoi(x[2])-strtoi(x[1])  ) ) 
  peaks_sizes = unlist(lapply(ends_all, FUN = function(x) strtoi(x[2])-strtoi(x[1])  )   )
                 
  # plot sizes of peaks
  hist(   small_peaks_sizes , breaks=500  )
  hist( peaks_sizes , breaks=1000 , xlim = c(0,5000) )
  abline(v=100)
  
  FRiP <- meta$atac_peak_region_fragments / meta$atac_fragments * 100
  hist(FRiP, breaks=500)
  #abline(v=3)
  abline(v=0.5)
  hist(FRiP, breaks=500, ylim = c(0,100))
  hist(FRiP[meta$is_cell==1] , breaks=500, ylim = c(0,100))
  
  eD_out = read.table(paste0("/mnt/beegfs6/home3/ahringer/em613/analysis/multiomics/emptryDrops_multiome2023_eDv3/data/output/realdata/", old_date,"/", stub, "/", stub, "_eD_multiome.tsv"), sep="\t", header=T)  
  named_FRiP <- FRiP
  names(named_FRiP) <- meta$barcode
  hist(named_FRiP[ eD_out$Row.names[  eD_out$FDR_multi < 0.001 & !is.na(eD_out$FDR_multi)]    ] , breaks=500, ylim = c(0,300))
  
  # plot eD and cR together
  c1 <- rgb(173,216,230, max = 255, alpha = 100, names = "darkblue")
  c2 <- rgb(255,192,203, max = 255, alpha = 100, names = "darkred")
  A = named_FRiP[ eD_out$Row.names[  eD_out$FDR_multi < 0.001 & !is.na(eD_out$FDR_multi)]    ]
  B = FRiP[meta$is_cell==1]
  b <- - 0.002 
  e <- max(c(A,B) ) + 5
  ax <- pretty(b:e, n = 300) 
  hg_eD =   hist(A , breaks=ax, plot = FALSE, border = F)
  hg_cR =   hist(B , breaks=ax, plot = FALSE, border=F)
  plot(hg_cR, col = c1, ylim=c(0,500), border = F) # Plot 1st histogram using a transparent color
  plot(hg_eD, col = c2, add = TRUE, ylim=c(0,500), border = F)
  legend(20, 500, lwd=3, col=c(c1, c2) , lty=1,
         legend=c( "cellRanger-arc", "EmptyDrops_multiome") )
  
  blue_trans = rgb(0, 0, 255, max = 255, alpha = 125, names = "blue50")
  hist(meta$atac_fragments[which(FRiP==100)], breaks = seq(0, max(meta$atac_fragments[which(FRiP==100)]), length.out = 201), xlim = c(0,100))
  abline(v=1, col=blue_trans)
  hist(meta$atac_fragments[which(FRiP==50)], breaks = seq(0, max(meta$atac_fragments[which(FRiP==50)]), length.out = 201), xlim = c(0,100))
  abline(v=2, col=blue_trans)
  hist(meta$atac_fragments[which(FRiP==25)], breaks = seq(0, max(meta$atac_fragments[which(FRiP==25)]), length.out = 301), xlim = c(0,50))
  abline(v=4, col=blue_trans)
  hist(meta$atac_fragments[which(FRiP==20)], breaks = seq(0, max(meta$atac_fragments[which(FRiP==20)]), length.out = 1001), xlim = c(0,100))
  abline(v=5, col=blue_trans)
  
  
  # PLOT OUTLINE FRiP after
  limits <- range(c(A,B))
  breaks <- seq(limits[1], limits[2], length.out=200)
  hg_eD =   hist(A , breaks=breaks, plot = FALSE)
  hg_cR =   hist(B , breaks=breaks, plot = FALSE)
  xrange <- range(hg_eD$mids)
  yrange <- range(hg_eD$counts)
  plot(0,0,type="n", xlim=xrange, ylim=yrange, xlab=expression("FRiP"), 
       ylab="Number of cells", cex.axis=1.2, cex.lab=1.4, main=stub, cex.main=1.4)
  shift <- 0
  colors <- c(c1, c2)
  modes = c(hg_eD$mids) #, hg_cR$mids)
  plotHistogramOutline(breaks+shift, hg_cR$counts, col="grey", lwd=2)
  shift <- 0
  plotHistogramOutline(breaks+shift, hg_eD$counts, col="salmon", lwd=2)
  legend("topright", col=c("grey","salmon"), legend=c( "cellRanger-arc", "EmptyDrops_multiome"), lwd=2, cex=1.2)
  #  hist outline end
  
  
  # PLOT OUTLINE FRiP after improved plot
  limits <- range(c(A,B))
  breaks <- seq(limits[1], limits[2], length.out=80)
  hg_eD =   hist(A , breaks=breaks, plot = FALSE)
  hg_cR =   hist(B , breaks=breaks, plot = FALSE)
  xrange <- range(hg_eD$mids)
  yrange <- range(hg_eD$counts)
  p_FRiP_after <- plot(0,0,type="n", xlim=xrange, ylim=yrange, xlab="", 
       ylab="", cex.axis=1.2, cex.lab=1.4, main="", font=2, cex.main=1.4)
  mtext(side=1, line=2, "FRiP", col="black", font=2,cex=1.5)
  mtext(side=2, line=2, "# of cells", col="black", font=2,cex=1.5)
  mtext(side=3, line=0.5, stub, col="black", font=2, cex=1.5)
  shift <- 0
  colors <- c(c1, c2)
  modes = c(hg_eD$mids) #, hg_cR$mids)
  plotHistogramOutline(breaks+shift, hg_cR$counts, col="grey", lwd=2)
  shift <- 0
  plotHistogramOutline(breaks+shift, hg_eD$counts, col="salmon", lwd=2)
  legend("topright", col=c("grey","salmon"), legend=c( "cellRanger-arc", "EmptyDrops_multiome"), lwd=2, cex=1.2)
  #  hist outline end  
  
  
  
  max_frip_of_excluded = max(FRiP[meta$excluded_reason==2]  )
  min(FRiP[meta$excluded_reason==2]  )
  min_frip_of_cells = min(FRiP[meta$is_cell==1]  )
  
  
  # FRiP before colored background
  hist(FRiP, ylim= c(0,10000), breaks=500)
  rect(-0.5, 0, max_frip_of_excluded, par("usr")[4], col = "pink", border = NA)
  par(new = TRUE)
  rect(max_frip_of_excluded, 0 ,100.5, par("usr")[4], col = "lightblue", border = NA)
  par(new = TRUE)
  hist(FRiP, ylim= c(0,10000), breaks=500)
  abline(v=max_frip_of_excluded, col="blue")
  abline(v=min_frip_of_cells, col="red")
  text(x=4,y=10000, round(max_frip_of_excluded,2) )
  
  # infer the padding used by cellranger to do cell calling
  eff_fgip = (max_frip_of_excluded + min_frip_of_cells) / 200
  eff_bpip = eff_fgip * 3*10^9
  eff_pad = (eff_bpip - sum(peaks_sizes)) / length(peaks_sizes)
  
  
  hist(FRiP[meta$excluded_reason==2] , breaks=500, ylim = c(0,1000))
  
  counts_over_fragments = colSums(count_matrix_w_nonpeak[seq(1, dim(count_matrix_w_nonpeak)[1]-1 ),]) / colSums(count_matrix_w_nonpeak)
  hist(counts_over_fragments, ylim= c(0,5000), breaks=500)
  hist(counts_over_fragments[meta$barcode[meta$is_cell==1]], breaks=500)
  #hist(counts_over_fragments[eD_cells], breaks=500, ylim=c(0,80))
  
  
  nCount = colSums(count_matrix_atac)
  max_counts_over_fragments_of_excluded = max( counts_over_fragments[ meta$barcode[meta$excluded_reason==2] ]  )
  min_counts_over_fragments_of_cells = min( counts_over_fragments[ meta$barcode[meta$is_cell==1] ]  )
  
  # By removing the small peaks how many more cells would be include?
  fr_reads_in_small_peaks = sum(small_peaks) / sum(meta$atac_fragments[meta$is_cell==1])
  fr_bp_in_small_peaks = ( sum(small_peaks_sizes)+eff_pad*length(small_peaks_sizes) ) / (3*10^9)
  max_num_more_cells <- sum( FRiP< eff_fgip*100 & FRiP > eff_fgip *100 - fr_bp_in_small_peaks*100 & !is.na(FRiP))
  FRiP_cR = FRiP[meta$is_cell==1]
  max_num_less_cR_cells <- sum(  FRiP_cR < eff_fgip *100 + fr_bp_in_small_peaks*100 & !is.na(FRiP_cR))
  plot.new() 
  text(x=0.5, y=.1, paste0("fr_reads_in_small_peaks=", fr_reads_in_small_peaks))
  text(x=0.5, y=.2, paste0("eff_pad=", eff_pad))
  text(x=0.5, y=0.3, paste0("fr_bp_in_small_peaks=", fr_bp_in_small_peaks))
  text(x=0.5, y=0.4, paste0("max_num_more_cells=", max_num_more_cells))
  text(x=0.5, y=0.5, paste0("max_num_less_cR_cells=", max_num_less_cR_cells))
  text(x=0.5, y=0.6, paste0("num_cR_cells=", sum(meta$is_cell)))

  hist(FRiP[meta$is_cell==1] , breaks=500, ylim = c(0,100))
  abline(v=eff_fgip *100 + fr_bp_in_small_peaks*100)
  hist(FRiP, breaks=500, ylim = c(0,10000))
  abline(v=max_frip_of_excluded, col="blue")
  abline(v=min_frip_of_cells, col="red")  

  dev.off()
}




