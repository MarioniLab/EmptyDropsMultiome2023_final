#setwd("~/analysis/EmptyDropsMultiome2023_final")
library(EmptyDropsMultiome)
library(DropletUtils)
library(Matrix)
library(rstudioapi)
library(Seurat)
library(mixtools)
library(tidyverse)
source("simulations/fcn_for_sim.R")
# library("optparse")
library(Seurat)
#library(Signac)
library(stringr)

current_date= paste(unlist(strsplit(as.character(Sys.Date()), "-")), collapse="")
opath <- paste0("data/output/figures/", current_date)
old_date = "20230711"

dir.create(opath,recursive=TRUE)

set.seed(73953024)

resdir <- paste0("data/output/sim/", old_date, "/results-sim")

eD_METADATA <- c(
  "pics-sim/pbmc_gran_sorted_eD_multiome_5000_5000-metadata.csv",
  "pics-sim/pbmc_gran_sorted_eD_multiome_5000_2000-metadata.csv",
  "pics-sim/pbmc_gran_sorted_eD_multiome_2000_5000-metadata.csv",
  "pics-sim/pbmc_gran_sorted_eD_multiome_2000_2000-metadata.csv"
)

eD_out <- c(
  "pics-sim/pbmc_gran_sorted_eD_multiome_5000_5000-output.csv",
  "pics-sim/pbmc_gran_sorted_eD_multiome_5000_2000-output.csv",
  "pics-sim/pbmc_gran_sorted_eD_multiome_2000_5000-output.csv",
  "pics-sim/pbmc_gran_sorted_eD_multiome_2000_2000-output.csv"
)


text_size =10



# figure 2e
f = file.path(resdir, "pbmc_gran_sortedeD-multiome-res.tsv")
tab <- read.table(f, header=TRUE, sep="\t", stringsAsFactors=FALSE)
tab$Method[tab$Method=="eD-multiome"]="eD"
tab$Method[tab$Method=="CellRanger-arc"]="cR"
tab$xcoords = seq_along(1:dim(tab)[1]) + cumsum(as.integer(tab$Method=="eD-multiome" & tab$G1Size==tab$G1Size[1] ))
observations2e = data.frame( "popul"= rep( c(rep("large cells", 4), rep("small cells", 4) ) , 4), 
                             "Recall"= c( tab$G1[1:4], tab$G2[1:4], 
                                          tab$G1[5:8], tab$G2[5:8],
                                          tab$G1[9:12], tab$G2[9:12],
                                          tab$G1[13:16], tab$G2[13:16]) ,
                             "xcoords"= c( seq(2:9), 9+seq(2:9), 18+seq(2:9), 27+seq(2:9)  ),
                             "method"= c( tab$Method[1:4], tab$Method[1:4], 
                                          tab$Method[5:8], tab$Method[5:8],
                                          tab$Method[9:12], tab$Method[9:12],
                                          tab$Method[13:16], tab$Method[13:16]) 
                             )


fig2e <- ggplot(observations2e, aes(x=xcoords, y=Recall, color=method )) + 
  geom_rect(color = NA, aes(xmin = 0, xmax = 8.5, ymin = 0, ymax = 1.01), fill = "lightgrey") + 
  geom_rect(color = NA, aes(xmin = 9.5, xmax = 17.5, ymin = 0, ymax = 1.01), fill = "lightgrey") + 
  geom_rect(color = NA, aes(xmin = 18.5, xmax = 26.5, ymin = 0, ymax = 1.01), fill = "lightgrey") + 
  geom_rect(color = NA, aes(xmin = 27.5, xmax = 36, ymin = 0, ymax = 1.01), fill = "lightgrey") + 
  geom_point(size=2, aes(shape=popul))+
  scale_color_manual(values=c("blue",  "green", "salmon", "orange") )+
  theme_bw()+
  xlab("")+
  theme(
    # Hide panel borders and remove grid lines
    panel.border = element_blank(),
    panel.background = element_rect("white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Change axis line
    axis.line = element_line(colour = "black"),
    axis.text=element_text(size=text_size),
    axis.title=element_text(size=text_size),
    #axis.text.x = element_blank(),
    #axis.text.x = c("2000", "2000"  ,"2000","2000"),
    legend.text = element_text(size=text_size),
    legend.title = element_blank(),
    axis.text.x = element_text(size=10, angle=30, vjust=.8, hjust=0.8)
  )  +
  scale_x_continuous(breaks=c(4, 13.5, 22.5, 31.5), labels= c("2000/2000", "2000/5000", "5000/2000", "5000/5000"))+
  guides(colour = guide_legend(override.aes = list(size=5)))






for (i in seq_along(eD_out) ) { 
  
  feDmeta = eD_METADATA[i]
  feDout = eD_out[i]
  
  eD_meta <- read.table( file.path("data/output/sim", old_date, feDmeta) , header = T, sep=",")
  eD.out_multi <- read.table( file.path("data/output/sim", old_date, feDout) , header = T, sep=",")
  
  stub <- substring(feDout, 30,47)
  g1 <- as.double( substring(stub, 10, 13) )
  g2 <- as.double( substring(stub, 15, 18) )
  
  ffile <- file.path(opath, paste0(stub, "_fig2.pdf"))  

  pdf(ffile, width=7.2, height=6.7)
  
  # fig 2a 
  eD.out_multi$identity_gt = "empty"
  eD.out_multi$identity_gt[substr(eD.out_multi$Row.names, 1,2)=="g2" ] = "small cells"
  eD.out_multi$identity_gt[substr(eD.out_multi$Row.names, 1,2)=="g1" ] = "large cells"
  
  observations2a = data.frame("log_atac"=log10(eD.out_multi$Total_chromatin+0.1),
                              "log_rna"=log10(eD.out_multi$Total_RNA+0.1),
                              "identity"=eD.out_multi$identity_gt   )
  
  fig2a <- ggplot(observations2a, aes(x=log_atac, y=log_rna, color=identity)) + 
    geom_point(size=0.1)+ 
    theme(
      # Hide panel borders and remove grid lines
      panel.border = element_blank(),
      panel.background = element_rect("white"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      # Change axis line
      axis.line = element_line(colour = "black"),
      legend.title = element_blank(),
      axis.text=element_text(size=text_size),
      axis.title=element_text(size=text_size)    #,legend.key.size=unit(10, 'cm')
    #) + scale_color_manual(values=c("deepskyblue2", "red", "purple") )+
    ) + scale_color_manual(values=c("cornflowerblue", "darkseagreen", "darkgoldenrod") )+
    ggplot2::ggtitle("ground truth of simulation")+
    guides(colour = guide_legend(override.aes = list(size=5)))
  
  
  # figure 2b
  eD.out_multi$identity[eD.out_multi$k_means==1 ] = "pass"
  eD.out_multi$identity[eD.out_multi$k_means==0 ] = "fail"
  observations2b = data.frame("log_atac"=log10(eD.out_multi$Total_chromatin+0.1),
                              "log_rna"=log10(eD.out_multi$Total_RNA+0.1),
                              "identity"=eD.out_multi$identity   )
  
  fig2b <- ggplot(observations2b, aes(x=log_atac, y=log_rna, color=identity)) + 
    geom_point(size=0.1)+ 
    theme(
      # Hide panel borders and remove grid lines
      panel.border = element_blank(),
      panel.background = element_rect("white"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.title = element_blank(),
      # Change axis line
      axis.line = element_line(colour = "black"),
      axis.text=element_text(size=text_size),
      axis.title=element_text(size=text_size)
    ) + scale_color_manual(values=c("red", "deepskyblue2") )+
    ggplot2::ggtitle("EmptyDropsMultiome on simulated data")+
    ggplot2::geom_abline(intercept = eD_meta$k_means_intercept,
                         slope = eD_meta$k_means_slope,
                        # linetype="dotted",
                        color = "black",
                        size=0.5)+
    guides(colour = guide_legend(override.aes = list(size=5)))

  
  # figure 2c
  eD.out_multi$identity = "fail"
  eD.out_multi$identity[ eD.out_multi$FDR_multi <= 0.001 & !is.na(eD.out_multi$FDR_multi) ] = "pass"
  observations2c = data.frame("log_atac"=log10(eD.out_multi$Total_chromatin+0.1),
                              "log_rna"=log10(eD.out_multi$Total_RNA+0.1),
                              "identity"=eD.out_multi$identity   )
  lines  = data.frame( name=c("higher", "k-means", "lower"),
                      slope=c(eD_meta$k_means_slope, eD_meta$k_means_slope, eD_meta$k_means_slope ),
                       intercept = c( eD_meta$higher_intercept, eD_meta$k_means_intercept,  eD_meta$lower_intercept ),
                      type=c("dashed", "solid", "dotted")
            )
  
  fig2c <- ggplot(observations2c, aes(x=log_atac, y=log_rna, color=identity)) + 
    geom_point(size=0.1)+ 
    theme(
      # Hide panel borders and remove grid lines
      panel.border = element_blank(),
      panel.background = element_rect("white"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      # Change axis line
      axis.line = element_line(colour = "black"),
      axis.text=element_text(size=text_size),
      legend.title = element_blank(),
      axis.title=element_text(size=text_size)
    ) + scale_color_manual(values=c("red", "deepskyblue2") )+
    ggplot2::ggtitle("EmptyDropsMultiome on simulated data")+
    guides(colour = guide_legend(override.aes = list(size=5)))
  
  fig2c <- fig2c + geom_abline(data = lines, aes(intercept=intercept, slope=slope, linetype=name))
    
    
  # figure 2d
  hard_cut <- log10(eD.out_multi$Total_RNA + 1) - eD_meta$k_means_slope * log10(eD.out_multi$Total_chromatin+1)
  emp.roc <- createRocPts(log(eD.out_multi$FDR_multi+1e-3), eD.out_multi$identity_gt)
  emp_rna.roc <- createRocPts(log(eD.out_multi$FDR_RNA+1e-3), eD.out_multi$identity_gt)
  lib.roc <- createRocPts(-hard_cut, eD.out_multi$identity_gt)

  emp.fdr = emp.roc[,"empty"]/rowSums(emp.roc)
  emp.tp1 = emp.roc[,"large cells"]/g1
  emp.tp2 = emp.roc[,"small cells"]/g2
  lib.fdr = lib.roc[,"empty"]/rowSums(lib.roc)
  lib.tp1 = lib.roc[,"large cells"]/g1
  lib.tp2 = lib.roc[,"small cells"]/g2
  emp_rna.fdr = emp_rna.roc[,"empty"]/rowSums(emp_rna.roc)
  emp_rna.tp1 = emp_rna.roc[,"large cells"]/g1
  emp_rna.tp2 = emp_rna.roc[,"small cells"]/g2

  length_rep = length(emp.fdr)+1
  observations2d = data.frame(   
    "FDR" = c( c(0, emp.fdr), c(0, emp.fdr), c(0, lib.fdr), c(0, lib.fdr), c(0, emp_rna.fdr), c(0, emp_rna.fdr) ),
    "TPR" = c( c(0, emp.tp1), c(0, emp.tp2), c(0, lib.tp1), c(0, lib.tp2), c(0, emp_rna.tp1), c(0, emp_rna.tp2) ),
    "identity" = c( rep("eD (large cells) ", length_rep), rep("eD (small cells)", length_rep), 
                    rep("cR (large cells)", length_rep), rep("cR (small cells)", length_rep), 
                    rep("eD-rna (large cells) ", length_rep), rep("eD-rna (small cells)", length_rep)     )
    )
  
  fig2d <- ggplot(observations2d, aes(x=FDR, y=TPR, colour=identity, lty = identity)) + 
    geom_line()+
    theme_bw()+
    scale_x_continuous(limits = c(0, 0.05))+
    scale_color_manual(values=c("blue", "blue", "salmon", "salmon", "orange", "orange") )+
    scale_linetype_manual(values = c("solid", "dashed", "solid", "dashed", "solid", "dashed") ) +
    guides(colour = guide_legend(override.aes = list(size=5)))+ 
    theme(
      # Hide panel borders and remove grid lines
      panel.border = element_blank(),
      panel.background = element_rect("white"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.title = element_blank(),
      # Change axis line
      axis.line = element_line(colour = "black"),
      axis.text=element_text(size=text_size),
      axis.title=element_text(size=text_size),
      axis.text.x = element_text(size=10, angle=30, vjust=.8, hjust=0.8)
    ) #+
    #scale_x_continuous(breaks=c(0, 0.01, 0.02), labels= c("0.00", "0.01", "0.02"))
    
  

  library(ggpubr)
  figure <- ggarrange(fig2a, fig2c, fig2e, fig2d, 
                      labels = c("A", "B", "C", "D"),
                      ncol = 2, nrow =2  )
  print(figure)
  
  
  
  
  dev.off()
  ffile <- file.path(opath, paste0(stub, "_fig2.jpeg"))  
  ggsave(ffile, plot = figure, width=7.2, height=6.0, units="in")
  
}






