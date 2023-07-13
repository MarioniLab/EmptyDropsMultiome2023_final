#setwd("/mnt/beegfs6/home3/ahringer/em613/analysis/multiomics/emptryDrops_multiome2023_eDv3")
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
library(grid)


current_date= paste(unlist(strsplit(as.character(Sys.Date()), "-")), collapse="")
opath <- paste0("data/output/figures/", current_date)
old_date = "20230608"

text_size = 10

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





for (i in seq_along(ALLFILES) ) { 
  
  fname = ALLFILES[i]
  fmeta = METADATA[i]
  feDmeta = eD_METADATA[i]
  
  # sce <- Read10X_h5( file.path("data/input", fname) )
  meta <- read.table( file.path("data/input", fmeta) , header = T, sep=",")
  eD_meta <- read.table( file.path("data/output/realdata", old_date, feDmeta) , header = T, sep="\t")
  
  stub <- sub("/.*", "", fname, "_qc")
  
  ffile <- file.path(opath, paste0(stub, "_fig1.pdf"))  
  
  # count_matrix_atac <- sce[["Peaks"]]
  # count_matrix_rna <- sce[["Gene Expression"]]
  # 
  
  pdf(ffile, width=7.2, height=6.7)
  
  FRiP <- meta$atac_peak_region_fragments / meta$atac_fragments * 100
  
  eD.out_multi = read.table(paste0("data/output/realdata/", old_date,"/", stub, "/", stub, "_eD_multiome.tsv"), sep="\t", header=T)  
  named_FRiP <- FRiP
  names(named_FRiP) <- meta$barcode
  
  max_frip_of_excluded = max(FRiP[meta$excluded_reason==2]  )
  min_frip_of_cells = min(FRiP[meta$is_cell==1]  )
  
  A = named_FRiP[ eD.out_multi$Row.names[  eD.out_multi$FDR_multi < 0.001 & !is.na(eD.out_multi$FDR_multi)]    ]
  B = FRiP[meta$is_cell==1]
  
  # FIGURE3a: FRiP after
  limits <- range(c(A,B))
  breaks <- seq(limits[1], limits[2], length.out=80)
  hg_eD =   hist(A , breaks=breaks, plot = FALSE)
  hg_cR =   hist(B , breaks=breaks, plot = FALSE)
  x <- rep(breaks, each=2)
  y <- c(0, rep(hg_cR$counts, each=2), 0)
  z <- c(0, rep(hg_eD$counts, each=2), 0)
  df =data.frame("xcoord" = c(x,x) , "ycoord" = c(y,z), "method"= c(rep("cR", length(y)), rep("eD", length(y)))   )
  fig2a <- ggplot(data=df, aes(x=xcoord, y=ycoord, group=method)) +
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
    )+scale_color_manual(values=c("grey", "salmon" ))
  fig2a
  
  
  
  
  
  # FIGURE 1a: FRiP before
  observations = data.frame("FRiP"=FRiP)
  fig1a <- ggplot2::ggplot(observations) + 
    geom_rect(aes(xmin = -1, xmax = max_frip_of_excluded, ymin = 0, ymax = 9000), fill = "pink") + 
    geom_rect(aes(xmin = max_frip_of_excluded, xmax = 101, ymin = 0, ymax = 9000), fill = "deepskyblue2") + 
    ggplot2::geom_histogram(ggplot2::aes(x = FRiP), binwidth = 0.1, colour = "black", 
                            fill = "black") +
    ggplot2::scale_x_continuous( limit = c(-0.2, 100.2), oob = function(x, limits) x)+
    ggplot2::scale_y_continuous( limit = c(0, 8000), oob = function(x, limits) x)+
    ggplot2::ylab("Frequency") +
    ggplot2::xlab("FRiP(%)") +
    ggplot2::ggtitle("cR FRiP threshold")+
    ggplot2::geom_vline(xintercept = min_frip_of_cells,
                        # linetype="dotted",
                        color = "blue",
                        size=0.5)+
    ggplot2::geom_vline(xintercept = max_frip_of_excluded,
                        # linetype="dotted",
                        color = "red",
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
      axis.title=element_text(size=text_size)
    )
  
  
  # FIGURE 1c: log_comb before
  log_comb <- log10(eD.out_multi$Total_RNA+0.1) - eD_meta[["k_means_slope"]] * log10(eD.out_multi$Total_chromatin+0.1)
  observations1b = data.frame("log_comb"=log_comb)
  slop = round(abs(eD_meta[["k_means_slope"]]), 2)
  fig1c <- ggplot2::ggplot(observations1b) + ggplot2::theme_bw()+
    ggplot2::theme_bw()+
    geom_rect(aes(xmin = -1, xmax = eD_meta[["k_means_intercept"]], ymin = 0, ymax = 9000), fill = "pink") + 
    geom_rect(aes(xmin = eD_meta[["k_means_intercept"]], xmax = 10, ymin = 0, ymax = 9000), fill = "deepskyblue2") + 
    ggplot2::geom_histogram(ggplot2::aes(x = log_comb), binwidth = 0.1, colour = "black", 
                            fill = "black") +
    ggplot2::scale_x_continuous( limit = c(-1, 10), oob = function(x, limits) x)+
    ggplot2::scale_y_continuous( limit = c(0, 8000), oob = function(x, limits) x)+
    ggplot2::ylab("Frequency") +
    #ggplot2::xlab(expression(Log[10]~"(ATAC count)"-"slope * "~Log[10]~"(RNA count)")) +
    #ggplot2::xlab(expression(Log[10]~"(ATAC)"+abs(eD_meta[["k_means_slope"]])+"*"~Log[10]~"(RNA)")) +
    ggplot2::xlab(bquote( "log(ATAC)+"~.(slop)~"* log(RNA)") ) +
    ggplot2::ggtitle("cR count threshold in 1d")+
    ggplot2::geom_vline(xintercept = eD_meta[["k_means_intercept"]],
                        # linetype="dotted",
                        color = "red",
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
      axis.title=element_text(size=text_size)
    )  
  
  
  # FIGURE 1b: k-means clustering before
  log_comb <- log10(eD.out_multi$Total_RNA) - eD_meta[["k_means_slope"]] * log10(eD.out_multi$Total_chromatin)
  observations1c = data.frame("log_atac"=log10(eD.out_multi$Total_chromatin+0.1),
                               "log_rna"=log10(eD.out_multi$Total_RNA+0.1),
                               "accepted"=log10(eD.out_multi$Total_chromatin+0.1)*eD_meta[["k_means_slope"]] + eD_meta[["k_means_intercept"]] > log10(eD.out_multi$Total_RNA+0.1)   )
  observations1c$accepted[observations1c$accepted=="FALSE"]="pass"
  observations1c$accepted[observations1c$accepted=="TRUE"]="rejected"
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
      axis.title=element_text(size=text_size),
      legend.title = element_blank()
    ) +
    scale_color_manual(values=c("deepskyblue2", "pink" ))+
    ggplot2::ggtitle("cR count threshold")+
    ggplot2::ylab("log(RNA)") +
    ggplot2::xlab("log(ATAC)")
  
  
  # FIGURE 1d: ATAC counts
  observations1d = data.frame("nCount"=eD.out_multi$Total_chromatin)
  vline2 = eD_meta[["lower_atac"]]
  vline1 = eD_meta[["barhop_atac"]]
  fig1d <- 
    ggplot2::ggplot(observations1d) + ggplot2::theme_bw() +
    ggplot2::geom_histogram(ggplot2::aes(x = nCount), binwidth = 1, colour = "black", 
                            fill = "black") +
    ggplot2::scale_x_continuous( limit = c(0, vline2+100), oob = function(x, limits) x)+
    ggplot2::scale_y_continuous( limit = c(0,50000 ), oob = function(x, limits) x)+
    ggplot2::ylab("Frequency") +
    ggplot2::xlab("ATAC Count") +
    ggplot2::ggtitle("Final GMM Fit")+
    # ggplot2::geom_vline(xintercept = vline1,
    #                     # linetype="dotted",
    #                     color = "purple",
    #                     size=0.5)+
    # ggplot2::geom_vline(xintercept = vline2,
    #                     # linetype="dotted",
    #                     color = "blue",
    #                     size=0.5)+ 
    theme(
      # Hide panel borders and remove grid lines
      panel.border = element_blank(),
      panel.background = element_rect("white"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      # Change axis line
      axis.line = element_line(colour = "black"),
      axis.text=element_text(size=text_size),
      #axis.title=element_text(size=text_size,face="bold")
      axis.title=element_text(size=text_size)
    )   
  
  
  # FIGURE 1e: RNA counts
  observations1e = data.frame("nCount"=eD.out_multi$Total_RNA)
  vline2 = eD_meta[["lower_rna"]]
  vline1 = eD_meta[["barhop_rna"]]
  fig1e <- 
    ggplot2::ggplot(observations1e) + ggplot2::theme_bw() +
    ggplot2::geom_histogram(ggplot2::aes(x = nCount), binwidth = 1, colour = "black", 
                            fill = "black") +
    ggplot2::scale_x_continuous( limit = c(0, vline2+100), oob = function(x, limits) x)+
    ggplot2::scale_y_continuous( limit = c(0,50000 ), oob = function(x, limits) x)+
    ggplot2::ylab("Frequency") +
    ggplot2::xlab("RNA Count") +
    ggplot2::ggtitle("Final GMM Fit")+
    # ggplot2::geom_vline(xintercept = vline1,
    #                     # linetype="dotted",
    #                     color = "purple",
    #                     size=0.5)+
    # ggplot2::geom_vline(xintercept = vline2,
    #                     # linetype="dotted",
    #                     color = "blue",
    #                     size=0.5) + 
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
  
  
  library(ggpubr)
  figure <- ggarrange(fig1a, fig1b, fig1c, fig1e, fig1d, 
                      labels = c("A", "B", "C", "D", "E"),
                      ncol = 2, nrow = 3, common.legend = TRUE, legend="bottom")
  print(figure)
  dev.off()

#   # pretty pdf version
#   library(grid)
#   # Move to a new page
#   grid.newpage()
#   # Create layout : nrow = 3, ncol = 2
#   pushViewport(viewport(layout = grid.layout(nrow = 3, ncol = 2)))
#   # A helper function to define a region on the layout
#   define_region <- function(row, col){
#     viewport(layout.pos.row = row, layout.pos.col = col)
#   }                              
#   pdf(file.path(opath, paste0(stub, "_fig1_pretty.pdf")) )                        
#   figure <- ggarrange(fig1b, fig1c, fig1d, fig1e,
#                     labels = c("B", "C", "D", "E"),
#                     ncol = 2, nrow = 2, common.legend = TRUE, legend="top")
#   print(fig1a, vp = define_region(row = 1, col = 1:2))   # Span over two columns
#   print(figure, vp = define_region(row = 2:3, col = 1:2))
#   dev.off() 
                
                                
  library("gridExtra")
  pdf(file.path(opath, paste0(stub, "_fig1_pretty.pdf")) )                        
  #grid.arrange(fig1a, fig1b, fig1c, layout_matrix = matrix(c(1, 3, 2, 3), nrow = 2))
  grid.arrange(fig1a, fig1b, fig1c, fig1e, fig1d, layout_matrix = rbind(c(1, 1),
                                                          c(2, 3),
                                                          c(4, 5))
              )
  dev.off() 
                                
                                
  ffile <- file.path(opath, paste0(stub, "_fig1.jpeg"))  
  ggsave(ffile, plot = figure, width=7.2, height=6.7, units="in")
                                
  
  ffile <- file.path(opath, paste0(stub, "_fig1_pretty.jpeg"))  
  ggsave(ffile, plot =   grid.arrange(fig1a, fig1b, fig1c, fig1e, fig1d, labels=c("A", "B", "C", "D", "E"), layout_matrix = rbind(c(1, 1),
                                                          c(2, 3),
                                                          c(4, 5))
              ), 
         width=7.2, height=6.7, units="in")
  
}


