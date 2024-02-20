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
library(eulerr)
library(ggridges)
library(ggpubr)


current_date= paste(unlist(strsplit(as.character(Sys.Date()), "-")), collapse="")
opath <- paste0("data/output/figures/", current_date)
old_date = "20230712" #"20230608"
less_old_date = "downstream_manybad_20230713" #"downstream_manybad_20230608"

dir.create(opath,recursive=TRUE)

set.seed(73953024)

ALLFILES <- c(
  "valentina_8176/FCA_GND10288176_raw_feature_bc_matrix.h5",
  "valentina_8177/FCA_GND10288177_raw_feature_bc_matrix.h5"
)

text_size=7


genes_germcells_to_mesGATA4 = c( "DAZL", "UPK3B", "GATA4", "LHX9", "NR5A1", "WNT6", "IRX3", "FOXL2", "ARX")
genes_mesGATA2_to_neural = c( "TCF21", "PDGFRA", "DCN", "GATA2", "NR2F1", "PDGFRB")
#genes_mesGATA2_to_neural2 = c( "MYH11", "PTPRC", "CDH5", "PAX8", "EPCAM", "HBA1", "ASCL1")
genes_mesGATA2_to_neural2 = c( "MYH11", "PTPRC", "CDH5", "PAX8", "EPCAM", "HBA1")
genes_PGC_FGC_FGCmitotic_oogoniaSTRA8 = c("DAZL", "IFITM1", "NANOG", "NANOS3", "POU5F1", "DDX4", "MAEL", "ZGLP1", "STRA8")
genes_oogoniameiotic_preoocyte_oocyte_prespermatogonia = c( "MEIOC", "SYCP1", "FIGLA", "LHX8", "NOBOX", "GDF9", "ZP3", "DCC")


ffile <- file.path(opath, "fig4.pdf") 
  

  # sample 1
  i=1
  fname = ALLFILES[i]
  stub <- sub("/.*", "", fname)
  srat <- readRDS(  file.path("data/output/realdata", old_date, stub, less_old_date, "srat_vale_clean.rds")   ) 
  srat_upstream <- readRDS(  file.path("data/output/realdata", old_date, stub, paste0("srat_", stub, ".rds")  ) ) 
  
  srat$seurat_clusters_old <- srat$seurat_clusters
  srat <- RenameIdents(object = srat, `10` = "ooc", `5` = "pre-ooc-2", `6` = "pre-ooc-1", `2`="oog-m", `3`="oog-st", `8`="PGC", 
                       `9`="CoEp", `4`="supp1", `0`="supp2",  `1`="supp3", `11`="supp3",
                       `7`="Mes",  `13`="Imm", `12`="Endo")
  srat$seurat_clusters <- srat@active.ident
  
  
  fill_var <- factor(srat$comparison, levels=c('cellRanger-arc', 'both', 'Emptydrops-multiome'))
  venn_df_ggplot <- data.frame("cluster"=srat$seurat_clusters,
                               "comparison"=fill_var
                               )

  barplot1 <- ggplot(data=venn_df_ggplot, aes( fill=comparison)) +
    geom_bar(mapping = aes(x=cluster), position = "fill")+
    scale_fill_manual(values=c("blue", "grey","red" ) )+
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
      axis.title=element_text(size=text_size),
      axis.text.x = element_text(size=text_size, angle=45, vjust=.8, hjust=0.8)
    ) + 
    xlab("cluster")+
    ylab("composition")+
    ggplot2::ggtitle("")+
    guides(colour = guide_legend(override.aes = list(size=5)))
  barplot1
  
  umap_comparison1 <- DimPlot(srat, group.by = "comparison", sizes.highlight=0.1, reduction = "umap",  cols = c("grey", "blue", "red"))+
    theme(
      # Hide panel borders and remove grid lines
      panel.border = element_blank(),
      panel.background = element_rect("white"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.text=element_text(size=text_size),
      # Change axis line
      axis.line = element_line(colour = "black"),
      axis.text=element_text(size=text_size),
      legend.title = element_blank(),
      axis.title=element_text(size=text_size),
      axis.text.x = element_text(size=text_size, angle=45, vjust=.8, hjust=0.8)
    ) +
    ggplot2::ggtitle("")
  umap_comparison1
  
  # srat$sub.cluster_old <- srat$sub.cluster
  # srat <- RenameIdents(object = srat, `10` = "ooc", `6` = "pre-ooc-2", `5` = "pre-ooc-1", `2`="oog-m", `3`="oog-st", `8`="PGC", 
  #                      `9`="CoEp", `0`="supp1", `11`="supp2", `4`="supp2", `1`="supp3",
  #                      `7`="Mes", `13`="Imm", `12`="Endo")
  # srat$sub.cluster <- srat@active.ident
  # 
  
  umap1 <- DimPlot(srat, group.by = "seurat_clusters", sizes.highlight=0.1, label.size = text_size/3, label=T, reduction = "umap")+guides(colour = "none")+
    theme(
      # Hide panel borders and remove grid lines
      panel.border = element_blank(),
      panel.background = element_rect("white"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.text=element_text(size=text_size),
      # Change axis line
      axis.line = element_line(colour = "black"),
      axis.text=element_text(size=text_size),
      legend.title = element_blank(),
      axis.title=element_text(size=text_size),
      axis.text.x = element_text(size=text_size, angle=45, vjust=.8, hjust=0.8)
    )+
    ggplot2::ggtitle("")
  umap1

  
  germ_cells1 <- DotPlot(srat, dot.scale=3, features = c(genes_PGC_FGC_FGCmitotic_oogoniaSTRA8, genes_oogoniameiotic_preoocyte_oocyte_prespermatogonia), 
                         assay = "SCT")+
    theme(
      # Hide panel borders and remove grid lines
      panel.border = element_blank(),
      panel.background = element_rect("white"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.text=element_text(size=text_size),
      legend.title = element_text(size=text_size),
      # Change axis line
      axis.line = element_line(colour = "black"),
      axis.text=element_text(size=text_size),
      axis.title=element_text(size=text_size),
      axis.text.x = element_text(size=text_size, angle=45, vjust=.8, hjust=0.8)
    ) + 
    ylab("clusters")+
    ggplot2::ggtitle("")+
    guides(colour = guide_legend(override.aes = list(size=5)))
  
  supporting_somatic1 <- DotPlot(srat, dot.scale=3, features = c(genes_germcells_to_mesGATA4, genes_mesGATA2_to_neural, genes_mesGATA2_to_neural2), 
                                 assay = "SCT")+
    theme(
      # Hide panel borders and remove grid lines
      panel.border = element_blank(),
      panel.background = element_rect("white"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.text=element_text(size=text_size),
      legend.title = element_text(size=text_size),
      # Change axis line
      axis.line = element_line(colour = "black"),
      axis.text=element_text(size=text_size),
      axis.title=element_text(size=text_size),
      axis.text.x = element_text(size=text_size, angle=45, vjust=.8, hjust=0.8)
    ) + 
    ylab("clusters")+
    ggplot2::ggtitle("")+
    guides(colour = guide_legend(override.aes = list(size=5)))
  
  max_frip_of_excluded = max(srat$FRiP[srat$excluded_reason==2]  )
  df_FRiP <- data.frame("frip"=srat$FRiP, "cluster"=srat$seurat_clusters)
  ridgeplot1 <- ggplot(df_FRiP, aes(x = frip, y = cluster)) + 
    geom_rect(aes(xmin = -0.01, xmax = max_frip_of_excluded, ymin = -Inf, ymax = Inf), fill = "pink") + 
    geom_rect(aes(xmin = max_frip_of_excluded, xmax = 35, ymin = -Inf, ymax = Inf), fill = "deepskyblue2") + 
    geom_density_ridges(stat = "binline", bins = 80, scale = 0.95, draw_baseline = FALSE)+
    theme_ridges(grid = FALSE, center_axis_labels = TRUE)+
    geom_vline(xintercept = max_frip_of_excluded)
  
  
  # euler figures
  listInput1 <- list(eD = colnames(srat)[srat@meta.data[["comparison"]] %in% c("Emptydrops-multiome", "both")]   , 
                    cR = colnames(srat)[srat@meta.data[["comparison"]] %in% c("cellRanger-arc", "both") ] 
  )
  overlaps1 <- euler(listInput1, shape = "ellipse")
  euler_figure1 <- plot(overlaps1, 
                       quantities = TRUE,
                       labels = list(font = 4))
  euler_figure1
  
  listInput1_bf_qc <- list(eD = colnames(srat_upstream)[srat_upstream@meta.data[["comparison"]] %in% c("Emptydrops-multiome", "both")]   , 
                     cR = colnames(srat_upstream)[srat_upstream@meta.data[["comparison"]] %in% c("cellRanger-arc", "both") ] 
  )
  overlaps1_bf_qc <- euler(listInput1_bf_qc, shape = "ellipse")
  euler_figure1_bf_qc <- plot(overlaps1_bf_qc, 
                        quantities = TRUE,
                        labels = list(font = 4))
  euler_figure1_bf_qc
  
  # quality statistics
  T_HQ_R_eD1 =  (overlaps1$original.values["eD"]+overlaps1$original.values["eD&cR"]) / sum(overlaps1$original.values)
  T_HQ_R_cR1 =  (overlaps1$original.values["cR"]+overlaps1$original.values["eD&cR"]) / sum(overlaps1$original.values)
  LQ_DR_eD1 = ( overlaps1_bf_qc$original.values["eD"]+overlaps1_bf_qc$original.values["eD&cR"]  - 
                  overlaps1$original.values["eD"]-overlaps1$original.values["eD&cR"]) / (   overlaps1$original.values["eD"]+overlaps1$original.values["eD&cR"]  )
  LQ_DR_cR1 = ( overlaps1_bf_qc$original.values["cR"]+overlaps1_bf_qc$original.values["eD&cR"]  - 
                  overlaps1$original.values["cR"]-overlaps1$original.values["eD&cR"]) /  (   overlaps1$original.values["cR"]+overlaps1$original.values["eD&cR"]  )
  quality_eD1 =  (   overlaps1$original.values["eD"]+overlaps1$original.values["eD&cR"]  ) /
    ( overlaps1_bf_qc$original.values["eD"]+overlaps1_bf_qc$original.values["eD&cR"]) 
  quality_cR1 =  (   overlaps1$original.values["cR"]+overlaps1$original.values["eD&cR"]  ) /
    ( overlaps1_bf_qc$original.values["cR"]+overlaps1_bf_qc$original.values["eD&cR"]) 
  print("the fraction of eD only droplets in oocytes is")
  print(sum(venn_df_ggplot$cluster=="ooc" & venn_df_ggplot$comparison=="Emptydrops-multiome") / sum(venn_df_ggplot$cluster=="ooc" ) )
  
  
  
  
  
  
  # -------------  second sample  ----------------------------------------------------------------------
  i=2
  fname = ALLFILES[i]
  stub <- sub("/.*", "", fname)
  srat <- readRDS(  file.path("data/output/realdata", old_date, stub, less_old_date, "srat_vale_clean.rds")   ) 
  srat_upstream <- readRDS(  file.path("data/output/realdata", old_date, stub, paste0("srat_", stub, ".rds")  ) ) 
  
  
  srat$seurat_clusters_old <- srat$seurat_clusters
  srat <- RenameIdents(object = srat, `10` = "ooc", `6` = "pre-ooc-2", `5` = "pre-ooc-1", `3`="oog-m", `4`="oog-st", `11`="oog-st", `8`="PGC", 
                       `9`="CoEp", `1`="supp1", `2`="supp2", `0`="supp3",
                       `7`="Mes",  `13`="Imm", `12`="Endo", `14`="Ery")
  srat$seurat_clusters <- srat@active.ident
  
  
  fill_var <- factor(srat$comparison, levels=c('cellRanger-arc', 'both', 'Emptydrops-multiome'))
  venn_df_ggplot <- data.frame("cluster"=srat$seurat_clusters,
                               "comparison"=fill_var
  )

  barplot2 <- ggplot(data=venn_df_ggplot, aes( fill=comparison)) +
    geom_bar(mapping = aes(x=cluster), position = "fill")+
    scale_fill_manual(values=c("blue", "grey","red" ) )+
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
      axis.title=element_text(size=text_size),
      axis.text.x = element_text(size=text_size, angle=45, vjust=.8, hjust=0.8)
    ) + 
    xlab("cluster")+
    ylab("composition")+
    ggplot2::ggtitle("")+
    guides(colour = guide_legend(override.aes = list(size=5)))
  barplot2
  
  umap_comparison2 <- DimPlot(srat, group.by = "comparison", sizes.highlight=0.1, reduction = "umap",  cols = c("grey", "blue", "red"))+
    theme(
      # Hide panel borders and remove grid lines
      panel.border = element_blank(),
      panel.background = element_rect("white"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.text=element_text(size=text_size),
      # Change axis line
      axis.line = element_line(colour = "black"),
      axis.text=element_text(size=text_size),
      legend.title = element_blank(),
      axis.title=element_text(size=text_size),
      axis.text.x = element_text(size=text_size, angle=45, vjust=.8, hjust=0.8)
    )+
    ggplot2::ggtitle("")
  umap_comparison2
  

  umap2 <- DimPlot(srat, group.by = "seurat_clusters", sizes.highlight=0.1, label.size = text_size/3, label=T, reduction = "umap")+guides(colour = "none")+
    theme(
      # Hide panel borders and remove grid lines
      panel.border = element_blank(),
      panel.background = element_rect("white"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.text=element_text(size=text_size),
      # Change axis line
      axis.line = element_line(colour = "black"),
      axis.text=element_text(size=text_size),
      legend.title = element_blank(),
      axis.title=element_text(size=text_size),
      axis.text.x = element_text(size=text_size, angle=45, vjust=.8, hjust=0.8)
    )+
    ggplot2::ggtitle("")
  umap2
  
  
  germ_cells2 <- DotPlot(srat, dot.scale=3, features = c(genes_PGC_FGC_FGCmitotic_oogoniaSTRA8, genes_oogoniameiotic_preoocyte_oocyte_prespermatogonia), 
                         assay = "SCT") +
    theme(
      # Hide panel borders and remove grid lines
      panel.border = element_blank(),
      panel.background = element_rect("white"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.text=element_text(size=text_size),
      legend.title = element_text(size=text_size),
      # Change axis line
      axis.line = element_line(colour = "black"),
      axis.text=element_text(size=text_size),
      axis.title=element_text(size=text_size),
      axis.text.x = element_text(size=text_size, angle=45, vjust=.8, hjust=0.8)
    ) + 
    ylab("clusters")+
    ggplot2::ggtitle("")+
    guides(colour = guide_legend(override.aes = list(size=5)))
  germ_cells2
  
  supporting_somatic2 <- DotPlot(srat, dot.scale=3, features = c(genes_germcells_to_mesGATA4, genes_mesGATA2_to_neural, genes_mesGATA2_to_neural2), 
                                 assay = "SCT")+
    theme(
      # Hide panel borders and remove grid lines
      panel.border = element_blank(),
      panel.background = element_rect("white"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.text=element_text(size=text_size),
      legend.title = element_text(size=text_size),
      # Change axis line
      axis.line = element_line(colour = "black"),
      axis.text=element_text(size=text_size),
      axis.title=element_text(size=text_size),
      axis.text.x = element_text(size=text_size, angle=45, vjust=.8, hjust=0.8)
    ) + 
    ylab("clusters")+
    ggplot2::ggtitle("")+
    guides(colour = guide_legend(override.aes = list(size=5)))
  
  max_frip_of_excluded = max(srat$FRiP[srat$excluded_reason==2]  )
  df_FRiP <- data.frame("frip"=srat$FRiP, "cluster"=srat$seurat_clusters)
  ridgeplot2 <-  ggplot(df_FRiP, aes(x = frip, y = cluster)) + 
    geom_rect(aes(xmin = -0.01, xmax = max_frip_of_excluded, ymin = -Inf, ymax = Inf), fill = "pink") + 
    geom_rect(aes(xmin = max_frip_of_excluded, xmax = 35, ymin = -Inf, ymax = Inf), fill = "deepskyblue2") + 
    geom_density_ridges(stat = "binline", bins = 80, scale = 0.95, draw_baseline = FALSE)+
    theme_ridges(grid = FALSE, center_axis_labels = TRUE)+
    geom_vline(xintercept = max_frip_of_excluded)
  
  # euler figures
  
  listInput2 <- list(eD = colnames(srat)[srat@meta.data[["comparison"]] %in% c("Emptydrops-multiome", "both")]   , 
                     cR = colnames(srat)[srat@meta.data[["comparison"]] %in% c("cellRanger-arc", "both") ] 
  )
  overlaps2 <- euler(listInput2, shape = "ellipse")
  euler_figure2 <- plot(overlaps2, 
                        quantities = TRUE,
                        labels = list(font = 4))
  euler_figure2
  
  listInput2_bf_qc <- list(eD = colnames(srat_upstream)[srat_upstream@meta.data[["comparison"]] %in% c("Emptydrops-multiome", "both")]   , 
                           cR = colnames(srat_upstream)[srat_upstream@meta.data[["comparison"]] %in% c("cellRanger-arc", "both") ] 
  )
  overlaps2_bf_qc <- euler(listInput2_bf_qc, shape = "ellipse")
  euler_figure2_bf_qc <- plot(overlaps2_bf_qc, 
                              quantities = TRUE,
                              labels = list(font = 4))
  euler_figure2_bf_qc
  
  # quality statistics
  T_HQ_R_eD2 =  (overlaps2$original.values["eD"]+overlaps2$original.values["eD&cR"]) / sum(overlaps2$original.values)
  T_HQ_R_cR2 =  (overlaps2$original.values["cR"]+overlaps2$original.values["eD&cR"]) / sum(overlaps2$original.values)
  LQ_DR_eD2 = ( overlaps2_bf_qc$original.values["eD"]+overlaps2_bf_qc$original.values["eD&cR"]  - 
                 overlaps2$original.values["eD"]-overlaps2$original.values["eD&cR"]) / (   overlaps2_bf_qc$original.values["eD"]+overlaps2_bf_qc$original.values["eD&cR"]  )
  LQ_DR_cR2 = ( overlaps2_bf_qc$original.values["cR"]+overlaps2_bf_qc$original.values["eD&cR"]  - 
                 overlaps2$original.values["cR"]-overlaps2$original.values["eD&cR"]) /  (   overlaps2_bf_qc$original.values["cR"]+overlaps2_bf_qc$original.values["eD&cR"]  )
  quality_eD2 =  (   overlaps2$original.values["eD"]+overlaps2$original.values["eD&cR"]  ) /
    ( overlaps2_bf_qc$original.values["eD"]+overlaps2_bf_qc$original.values["eD&cR"]) 
  quality_cR2 =  (   overlaps2$original.values["cR"]+overlaps2$original.values["eD&cR"]  ) /
    ( overlaps2_bf_qc$original.values["cR"]+overlaps2_bf_qc$original.values["eD&cR"]) 
  print("the fraction of eD only droplets in oocytes is")
  print(sum(venn_df_ggplot$cluster=="ooc" & venn_df_ggplot$comparison=="Emptydrops-multiome") / sum(venn_df_ggplot$cluster=="ooc" ) )
  
  
  
  
  
  
# figure 4  
pdf(ffile, width=7.2, height=6.7)
figure <- ggarrange(umap1, umap2, umap_comparison1, umap_comparison2, barplot1, barplot2,
                    labels = c("A", "B", "C", "D", "E", "F"),
                    ncol = 2, nrow = 3, common.legend = TRUE, legend="bottom", widths = c(1, 1, 1), heights = c(1,1,0.7)  )
print(figure)
dev.off()
ffile <- file.path(opath, "fig4.jpeg") 
ggsave(ffile, plot = figure, width=7.2, height=6.7, units="in")

ffile <- file.path(opath, "fig4A.pdf") 
ggsave(ffile, plot = umap1, width=7.2, height=6.7, units="in")
ffile <- file.path(opath, "fig4B.pdf") 
ggsave(ffile, plot = umap2, width=7.2, height=6.7, units="in")
ffile <- file.path(opath, "fig4C.pdf") 
ggsave(ffile, plot = umap_comparison1, width=7.2, height=6.7, units="in")
ffile <- file.path(opath, "fig4D.pdf") 
ggsave(ffile, plot = umap_comparison2, width=7.2, height=6.7, units="in")
ffile <- file.path(opath, "fig4E.pdf") 
ggsave(ffile, plot = barplot1, width=7.2, height=6.7, units="in")
ffile <- file.path(opath, "fig4F.pdf") 
ggsave(ffile, plot = barplot2, width=7.2, height=6.7, units="in")



# supp 1
ffile <- file.path(opath, "supp1.pdf") 
pdf(ffile, width=7.2, height=6.7)
figure_supp <- ggarrange(ridgeplot1, ridgeplot2,
                         labels = c("A", "B"),
                         ncol = 2, nrow = 1)
print(figure_supp)
dev.off()
ffile <- file.path(opath, "supp1.jpeg")  
ggsave(ffile, plot = figure_supp, width=7.2, height=6.7, units="in")

ffile <- file.path(opath, "supp1A.pdf")  
ggsave(ffile, plot = ridgeplot1, width=7.2, height=6.7, units="in")
ffile <- file.path(opath, "supp1B.pdf")  
ggsave(ffile, plot = ridgeplot2, width=7.2, height=6.7, units="in")



# supp 2
ffile <- file.path(opath, "supp2.pdf") 
pdf(ffile, width=7.2, height=6.7)
figure_supp <- ggarrange(umap1, umap2, germ_cells1, germ_cells2, supporting_somatic1, supporting_somatic2,
                    labels = c("A", "B", "C", "D", "E", "F"),
                    ncol = 2, nrow = 3, common.legend = TRUE, legend="bottom", widths = c(1, 1, 1), heights = c(0.7,1,1)  )
print(figure_supp)
dev.off()
ffile <- file.path(opath, "supp2.jpeg")  
ggsave(ffile, plot = figure_supp, width=7.2, height=6.7, units="in")

ffile <- file.path(opath, "supp2A.pdf")  
ggsave(ffile, plot = umap1, width=7.2, height=6.7, units="in")
ffile <- file.path(opath, "supp2B.pdf")  
ggsave(ffile, plot = umap2, width=7.2, height=6.7, units="in")
ffile <- file.path(opath, "supp2C.pdf")  
ggsave(ffile, plot = germ_cells1, width=7.2, height=6.7, units="in")
ffile <- file.path(opath, "supp2D.pdf")  
ggsave(ffile, plot = germ_cells2, width=7.2, height=6.7, units="in")
ffile <- file.path(opath, "supp2E.pdf")  
ggsave(ffile, plot = supporting_somatic1, width=7.2, height=6.7, units="in")
ffile <- file.path(opath, "supp2F.pdf")  
ggsave(ffile, plot = supporting_somatic2, width=7.2, height=6.7, units="in")



# supplemental figure: venns before and venns after
# supp 3
ffile <- file.path(opath, "supp3.pdf") 
pdf(ffile, width=7.2, height=6.7)
figure_supp <- ggarrange( euler_figure1_bf_qc, euler_figure2_bf_qc, euler_figure1, euler_figure2,
                         labels = c("Sample A before QC", "Sample B before QC", "Sample A after QC", "Sample B after QC"),
                         ncol = 2, nrow = 2 )
print(figure_supp)
dev.off()
ffile <- file.path(opath, "supp3.jpeg")  
ggsave(ffile, plot = figure_supp, width=7.2, height=6.7, units="in")

ffile <- file.path(opath, "supp3A.pdf")  
ggsave(ffile, plot = euler_figure1_bf_qc, width=7.2, height=6.7, units="in")
ffile <- file.path(opath, "supp3B.pdf")  
ggsave(ffile, plot = euler_figure2_bf_qc, width=7.2, height=6.7, units="in")
ffile <- file.path(opath, "supp3C.pdf")  
ggsave(ffile, plot = euler_figure1, width=7.2, height=6.7, units="in")
ffile <- file.path(opath, "supp3D.pdf")  
ggsave(ffile, plot = euler_figure2, width=7.2, height=6.7, units="in")










