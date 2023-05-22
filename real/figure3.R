setwd("/mnt/beegfs6/home3/ahringer/em613/analysis/multiomics/emptryDrops_multiome2023_eDv3")
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
library(eulerr)

current_date= paste(unlist(strsplit(as.character(Sys.Date()), "-")), collapse="")
opath <- paste0("data/output/figures/", current_date)
old_date = "20230329"
less_old_date = "downstream_manybad_20230413"

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
  "valentina_8177/FCA_GND10288177_raw_feature_bc_matrix.h5"
  
)


METADATA <- c(
  "valentina_8176/FCA_GND10288176_per_barcode_metrics.csv",
  "valentina_8177/FCA_GND10288177_per_barcode_metrics.csv"
)


eD_METADATA <- c(
  "valentina_8176/eD_metadata_valentina_8176.tsv",
  "valentina_8176/eD_metadata_valentina_8176.tsv"
)


text_size=10


for (i in seq_along(ALLFILES) ) { 
  
  fname = ALLFILES[i]
  fmeta = METADATA[i]
  feDmeta = eD_METADATA[i]
  
  meta <- read.table( file.path("data/input", fmeta) , header = T, sep=",")
  eD_meta <- read.table( file.path("data/output/realdata", old_date, feDmeta) , header = T, sep="\t")
  stub <- sub("/.*", "", fname)

  srat <- readRDS(  file.path("data/output/realdata", old_date, stub, less_old_date, "srat_vale_clean.rds")   ) 

  ffile <- file.path(opath, paste0(stub, "_fig3.pdf"))  
  
  pdf(ffile, width=7.2, height=6.7)
  
  
  FRiP <- meta$atac_peak_region_fragments / meta$atac_fragments * 100
  
  eD.out_multi = read.table(paste0("data/output/realdata/", old_date,"/", stub, "/", stub, "_eD_multiome.tsv"), sep="\t", header=T)  


  #A = named_FRiP[ eD.out_multi$Row.names[  eD.out_multi$FDR_multi < 0.001 & !is.na(eD.out_multi$FDR_multi)]    ]
  A = srat$FRiP[ srat@meta.data[["comparison"]] %in% c("Emptydrops-multiome", "both")    ]
  #B = FRiP[meta$is_cell==1]
  B = srat$FRiP[ srat@meta.data[["comparison"]] %in% c("cellRanger-arc", "both")  ]
  
  # FIGURE3a: FRiP after
  limits <- range(c(A,B))
  breaks <- seq(limits[1], limits[2], length.out=80)
  hg_eD =   hist(A , breaks=breaks, plot = FALSE)
  hg_cR =   hist(B , breaks=breaks, plot = FALSE)
  x <- rep(breaks, each=2)
  y <- c(0, rep(hg_cR$counts, each=2), 0)
  z <- c(0, rep(hg_eD$counts, each=2), 0)
  df =data.frame("xcoord" = c(x,x) , "ycoord" = c(y,z), "method"= c(rep("cR", length(y)), rep("eD", length(y)))   )
  frip_after <- ggplot(data=df, aes(x=xcoord, y=ycoord, group=method)) +
    geom_line( aes(color=method) )+ theme(
      # Hide panel borders and remove grid lines
      #panel.border = element_blank(),
      panel.background = element_rect("white"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      # Change axis line
      axis.line = element_line(colour = "black"),
      axis.text=element_text(size=text_size),
      axis.title=element_text(size=text_size)
    )+
    #scale_y_continuous(limit = c(1.8, 40), oob = function(x, limits) x)+
    scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) + 
    scale_y_continuous(expand = c(0, 0), limits = c(0, NA))+
    scale_color_manual(values=c("grey", "salmon" ))+
    ggplot2::ylab("Frequency") +
    ggplot2::xlab("FRiP") 
  frip_after

  
  #  cell_calling figure
  eD.out_multi$identity = "fail"
  eD.out_multi$identity[ eD.out_multi$Row.names %in% colnames(srat)[srat@meta.data[["comparison"]] %in% c("Emptydrops-multiome", "both")    ] ] = "pass"
  observations2c = data.frame("log_atac"=log10(eD.out_multi$Total_chromatin+0.1),
                              "log_rna"=log10(eD.out_multi$Total_RNA+0.1),
                              "identity"=eD.out_multi$identity   )
  
  cell_calling <- ggplot(observations2c, aes(x=log_atac, y=log_rna, color=identity)) + 
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
    ) + 
    xlab("log_atac")+
    ylab("log_rna")+
    scale_color_manual(values=c("red", "deepskyblue2") )+
    ggplot2::ggtitle("")+
    ggplot2::geom_abline(intercept = eD_meta$k_means_intercept,
                         slope = eD_meta$k_means_slope,
                         # linetype="dotted",
                         color = "black",
                         size=0.5)+
    ggplot2::geom_abline(intercept = eD_meta$lower_intercept,
                         slope = eD_meta$k_means_slope,
                         # linetype="dotted",
                         color = "black",
                         size=0.5)+
    ggplot2::geom_abline(intercept = eD_meta$higher_intercept,
                         slope = eD_meta$k_means_slope,
                         # linetype="dotted",
                         color = "black",
                         size=0.5)+
    guides(colour = guide_legend(override.aes = list(size=5)))
  cell_calling
  

  
  
  # FIGURE 3c: log_comb_after
  log_comb_after <- log10(eD.out_multi$Total_RNA+0.1) - eD_meta[["k_means_slope"]] * log10(eD.out_multi$Total_chromatin+0.1)
  names(log_comb_after) = eD.out_multi$Row.names
  slop = round(abs(eD_meta[["k_means_slope"]]), 2)
  A = log_comb_after[  colnames(srat)[srat@meta.data[["comparison"]] %in% c("Emptydrops-multiome", "both")    ]   ]
  B = log_comb_after[ colnames(srat)[srat@meta.data[["comparison"]] %in% c("cellRanger-arc", "both")    ]   ]
  
  limits <- range(c(A,B))
  breaks <- seq(limits[1], limits[2], length.out=60)
  hg_eD =   hist(A , breaks=breaks, plot = FALSE)
  hg_cR =   hist(B , breaks=breaks, plot = FALSE)
  x <- rep(breaks, each=2)
  y <- c(0, rep(hg_cR$counts, each=2), 0)
  z <- c(0, rep(hg_eD$counts, each=2), 0)
  df =data.frame("xcoord" = c(x,x) , "ycoord" = c(y,z), "method"= c(rep("cR", length(y)), rep("eD", length(y)))   )
  log_comb_after <- ggplot(data=df, aes(x=xcoord, y=ycoord, group=method)) +
    geom_line( aes(color=method) )+ theme(
      # Hide panel borders and remove grid lines
      #panel.border = element_blank(),
      panel.background = element_rect("white"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      # Change axis line
      axis.line = element_line(colour = "black"),
      axis.text=element_text(size=text_size),
      axis.title=element_text(size=text_size)
    )+
    scale_color_manual(values=c("grey", "salmon" ))+
    ggplot2::ylab("Frequency") +
    ggplot2::xlab(bquote( "log(ATAC)+"~.(slop)~"* log(RNA)") ) +
    ggplot2::ggtitle("cR count threshold in 1d")+    
    ggplot2::geom_vline(xintercept = eD_meta[["k_means_intercept"]],
                # linetype="dotted",
                color = "black",
                size=0.5)
  log_comb_after
  
  
  
  # FIGURE 1b: k-means clustering before
  log_comb <- log10(eD.out_multi$Total_RNA) - eD_meta[["k_means_slope"]] * log10(eD.out_multi$Total_chromatin)
  observations1c = data.frame("log_atac"=log10(eD.out_multi$Total_chromatin+0.1),
                               "log_rna"=log10(eD.out_multi$Total_RNA+0.1),
                               "accepted"=log10(eD.out_multi$Total_chromatin+0.1)*eD_meta[["k_means_slope"]] + eD_meta[["k_means_intercept"]] > log10(eD.out_multi$Total_RNA+0.1)   )
  fig1b <- ggplot(observations1c, aes(x=log_atac, y=log_rna, color=accepted)) + 
    geom_point()+ 
    theme(
      # Hide panel borders and remove grid lines
      panel.border = element_blank(),
      panel.background = element_rect("white"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      # Change axis line
      axis.line = element_line(colour = "black"),
      axis.text=element_text(size=text_size),
      axis.title=element_text(size=text_size)
    ) +
    scale_color_manual(values=c("deepskyblue2", "pink" ))+
    ggplot2::ggtitle("cR count threshold")
    
  
  
  # FIGURE 1d: ATAC counts
  observations3d = data.frame("nCount"=eD.out_multi$Total_chromatin)
  vline2 = eD_meta[["lower_atac"]]
  vline1 = eD_meta[["barhop_atac"]]
  ATAC_counts <- 
    ggplot2::ggplot(observations3d) + ggplot2::theme_bw() +
    ggplot2::geom_histogram(ggplot2::aes(x = nCount), binwidth = 1, colour = "black", 
                            fill = "black") +
    ggplot2::scale_x_continuous( limit = c(0, vline2+100), oob = function(x, limits) x)+
    ggplot2::scale_y_continuous(labels = scales::scientific, limit = c(0,50000 ), oob = function(x, limits) x)+
    ggplot2::ylab("Frequency") +
    ggplot2::xlab("ATAC Count") +
    ggplot2::ggtitle("Final GMM Fit")+
    ggplot2::geom_vline(xintercept = vline1,
                        # linetype="dotted",
                        color = "purple",
                        size=0.5)+
    ggplot2::geom_vline(xintercept = vline2,
                        # linetype="dotted",
                        color = "blue",
                        size=0.5)+ 
    theme(
      # Hide panel borders and remove grid lines
      panel.border = element_blank(),
      panel.background = element_rect("white"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      # Change axis line
      axis.line = element_line(colour = "black"),
      axis.text=element_text(size=text_size),
      #axis.title=element_text(size=15,face="bold")
      axis.title=element_text(size=text_size)
    )   
  
  
  # FIGURE 3e: RNA counts
  observations3e = data.frame("nCount"=eD.out_multi$Total_RNA)
  vline2 = eD_meta[["lower_rna"]]
  vline1 = eD_meta[["barhop_rna"]]
  RNA_counts <- 
    ggplot2::ggplot(observations3e) + ggplot2::theme_bw() +
    ggplot2::geom_histogram(ggplot2::aes(x = nCount), binwidth = 1, colour = "black", 
                            fill = "black") +
    ggplot2::scale_x_continuous( limit = c(0, vline2+100), oob = function(x, limits) x)+
    ggplot2::scale_y_continuous(labels = scales::scientific, limit = c(0,50000 ), oob = function(x, limits) x)+
    ggplot2::ylab("Frequency") +
    ggplot2::xlab("RNA Count") +
    ggplot2::ggtitle("Final GMM Fit")+
    ggplot2::geom_vline(xintercept = vline1,
                        # linetype="dotted",
                        color = "purple",
                        size=0.5)+
    ggplot2::geom_vline(xintercept = vline2,
                        # linetype="dotted",
                        color = "blue",
                        size=0.5) + 
    theme(
      # Hide panel borders and remove grid lines
      panel.border = element_blank(),
      panel.background = element_rect("white"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      # Change axis line
      axis.line = element_line(colour = "black"),
      axis.text=element_text(size=text_size),
      axis.title=element_text(size=text_size)
    )  
  RNA_counts
  
  
  
  # euler figure
  listInput <- list(eD = colnames(srat)[srat@meta.data[["comparison"]] %in% c("Emptydrops-multiome", "both")]   , 
                    cR = colnames(srat)[srat@meta.data[["comparison"]] %in% c("cellRanger-arc", "both") ] 
  )
  overlaps <- euler(listInput, shape = "ellipse")
  euler_figure <- plot(overlaps, 
             quantities = TRUE,
             labels = list(font = 4))
  euler_figure
  
  library(ggpubr)
  figure <- ggarrange(ATAC_counts, RNA_counts, euler_figure,  frip_after, cell_calling,  log_comb_after,
                      labels = c("A", "B", "C", "D", "E"),
                      ncol = 2, nrow = 3)
  print(figure)
  
  
  dev.off()
  
  ffile <- file.path(opath, paste0(stub, "_fig3.jpeg"))  
  ggsave(ffile, plot = figure, width=7.2, height=6.7, units="in")
  
}


