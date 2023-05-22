setwd("/mnt/beegfs6/home3/ahringer/em613/analysis/multiomics/emptryDrops_multiome2023_eDv3")
library(scDblFinder)
library(Seurat)
library(Signac)
library(scater)
library("optparse")
library(dplyr)
library(ggridges)
library(eulerr)
source("real/clustering_comparison.R")
source("simulations/fcn_for_sim.R")

# read date from user input
option_list = list(
  make_option(c("-d", "--date"), type="character", default=NULL,
              help="date of creation of seurat object aka location of directory", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
old_date <- "20230517"   #"20230329" #opt$date
samples <- c("valentina_8176", "valentina_8177")  #, "valentina_8177")

for (sample in samples){

# create directory
current_date= paste(unlist(strsplit(as.character(Sys.Date()), "-")), collapse="")
downpath <- paste0("data/output/realdata/", old_date, "/", sample, "/downstream_manybad_", current_date)
dir.create(downpath, recursive=TRUE)


fdoublet <- file.path(downpath, "downstream.pdf")
markers_tsv <- file.path(downpath, "markers.tsv")
markers_tsv_in_atac <- file.path(downpath, "markers_in_atac.tsv")
srat_vale_clean <- file.path(downpath, "srat_vale_clean.rds")
srat_atac_vale_clean <- file.path(downpath, "srat_atac_vale_clean.rds")

vale <- readRDS( paste0("data/output/realdata/", old_date, "/", sample, "/srat_", sample, ".rds") )
eD.out_multi <- read.table( paste0("data/output/realdata/", old_date, "/", sample,"/",sample, "_eD_multiome.tsv"), sep="\t", header=T) 

sce <- as.SingleCellExperiment(vale)
set.seed(123)
sce.mam.dbl <- scDblFinder(sce, clusters=colData(sce)$ident)
pdf(fdoublet)
print(plotUMAP(sce.mam.dbl, colour_by="scDblFinder.score"))
print(DimPlot(vale, label=T))
vale$scDblFinder.score  <- sce.mam.dbl$scDblFinder.score
t = table(sce.mam.dbl$scDblFinder.class)
print(FeaturePlot(vale, features=c("scDblFinder.score"),  order=T)+ggtitle(paste0("singlets: ", unname(t)[1], "; doublets: ", unname(t)[2]) ) )
print(VlnPlot(vale, features=c("scDblFinder.score"), y.max = 0.10)+stat_summary(fun.y = median, geom='point', size = 2, colour = "blue")+NoLegend()+ggtitle(paste0("singlets: ", unname(t)[1], "; doublets: ", unname(t)[2]) ) )

# remove clusters and doublets
if (sample=="valentina_8176"){
  #bad_clusters = c("3", "21")
  bad_clusters = c("3", "13", "20", "21", "22")
} else if (sample=="valentina_8177") {
  #bad_clusters = c("7", "17")
  #bad_clusters = c("7", "12", "17")
  bad_clusters = c("7", "11", "15")
}
keep_cells <- colnames(vale)[ !(vale$seurat_clusters %in% bad_clusters) & sce.mam.dbl$scDblFinder.class=="singlet" & vale$FRiP>median(vale$FRiP)-mad(vale$FRiP)   ] 
vale_clean <- subset(vale, cells = keep_cells)

plot.new()
rejected_cells = sum( vale$FRiP<median(vale$FRiP)-mad(vale$FRiP) )
text(x=0.2, y=.1, paste0("rejected_cells below 1 MAD less than median FRiP=", rejected_cells))


plot( vale$seurat_clusters[!vale$k_means], xlab="cluster", main="before qc: cells below k_means" )
plot( vale$seurat_clusters , xlab="cluster", main="before qc: total cells" )

# re-preprocess
vale_clean <- SCTransform(vale_clean)
vale_clean <- RunPCA(vale_clean, seed.use=42, features = VariableFeatures(object = vale_clean))
print(ElbowPlot(vale_clean, ndims = 50)     )
vale_clean <- FindNeighbors(vale_clean, dims = 1:50)
vale_clean <- FindClusters(vale_clean, resolution = 1, random.seed = 0)
vale_clean <- RunUMAP(vale_clean, dims = 1:50, seed.use = 42)

plot( vale_clean$seurat_clusters[!vale_clean$k_means], xlab="cluster", main="after qc: cells below k_means" )
plot( vale_clean$seurat_clusters , xlab="cluster", main="after qc: total cells" )


# subcluster the CoelEpith+Neuronal
vale_clean <- FindSubCluster(
  vale_clean,
  cluster=9,
  graph.name="SCT_snn",
  subcluster.name = "sub.cluster",
  resolution = 0.2,
  algorithm = 1
)
# vale_clean$sub.cluster[vale_clean$sub.cluster=="9_0"] = "9"
# vale_clean$sub.cluster[vale_clean$sub.cluster=="9_1"] = "14"
vale_clean$sub.cluster[vale_clean$sub.cluster=="9_0"] = "9"
vale_clean$sub.cluster[vale_clean$sub.cluster=="9_1"] = "9"


# DimPlot
print(DimPlot(vale_clean, label=T, group.by = "sub.cluster"))
print(VlnPlot(vale_clean,group.by = "sub.cluster", features="percent.mt")+ stat_summary(fun.y = median, geom='point', size = 10, colour = "blue") )
print(VlnPlot(vale_clean,group.by = "sub.cluster", features="percent.ribo")+ stat_summary(fun.y = median, geom='point', size = 10, colour = "blue"))
print(VlnPlot(vale_clean,group.by = "sub.cluster", features="nFeature_RNA")+ stat_summary(fun.y = median, geom='point', size = 10, colour = "blue"))
print(VlnPlot(vale_clean,group.by = "sub.cluster", features="nCount_RNA")+ stat_summary(fun.y = median, geom='point', size = 10, colour = "blue"))
print(VlnPlot(vale_clean,group.by = "sub.cluster", features="nCount_ATAC")+ stat_summary(fun.y = median, geom='point', size = 10, colour = "blue"))
print(VlnPlot(vale_clean,group.by = "sub.cluster", features="nFeature_ATAC")+ stat_summary(fun.y = median, geom='point', size = 10, colour = "blue"))

# FRiP plot per cluster
df_FRiP <- data.frame("frip"=vale_clean$FRiP, "cluster"=vale_clean$sub.cluster)
max_frip_of_excluded = max(vale_clean$FRiP[vale_clean$excluded_reason==2]  )
print(ggplot(df_FRiP, aes(x = frip, y = cluster, height = stat(density))) + 
        geom_density_ridges(stat = "binline", bins = 80, scale = 0.95, draw_baseline = FALSE)+
        theme_ridges(grid = FALSE, center_axis_labels = TRUE)+
        geom_vline(xintercept = max_frip_of_excluded))

# germ cells
genes_PGC_FGC_FGCmitotic_oogoniaSTRA8 = c("DAZL", "IFITM1", "NANOG", "NANOS3", "POU5F1", "DDX4", "MAEL", "ZGLP1", "STRA8")
genes_oogoniameiotic_preoocyte_oocyte_prespermatogonia = c( "MEIOC", "SYCP1", "FIGLA", "LHX8", "NOBOX", "ZP3", "GDF9", "SPOCD1", "MORC1")
print(DimPlot(vale_clean, group.by = "sub.cluster", label=T) | FeaturePlot(vale_clean, order=T, features = genes_PGC_FGC_FGCmitotic_oogoniaSTRA8 )& scale_colour_gradientn(colours = c( "cyan", "deepskyblue", "forestgreen", "darkorange", "darkred")) )
print(DimPlot(vale_clean, group.by = "sub.cluster", label=T) | FeaturePlot(vale_clean, order=T, features = genes_oogoniameiotic_preoocyte_oocyte_prespermatogonia )& scale_colour_gradientn(colours = c( "cyan", "deepskyblue", "forestgreen", "darkorange", "darkred")))
print(VlnPlot(vale_clean, group.by = "sub.cluster", features = genes_PGC_FGC_FGCmitotic_oogoniaSTRA8 ) )
print(VlnPlot(vale_clean, group.by = "sub.cluster", features = genes_oogoniameiotic_preoocyte_oocyte_prespermatogonia ))
print(VlnPlot(vale_clean, group.by = "sub.cluster", features = c("DCC", "ZP3") ))

# developing ovary
genes_germcells_to_mesGATA4 = c( "DAZL", "UPK3B", "GATA4", "LHX9", "NR5A1", "WNT6", "IRX3", "FOXL2", "ARX")
genes_mesGATA2_to_neural = c( "TCF21", "PDGFRA", "DCV", "GATA2", "NR2F1", "PDGFRB")
genes_mesGATA2_to_neural2 = c( "MYH11", "PTPRC", "CDH5", "PAX8", "EPCAM", "HBA1", "ASCL1")
print(DimPlot(vale_clean, group.by = "sub.cluster", label=T) | FeaturePlot(vale_clean, order=T, features = genes_germcells_to_mesGATA4 )& scale_colour_gradientn(colours = c( "cyan", "deepskyblue", "forestgreen", "darkorange", "darkred")) )
print(DimPlot(vale_clean, group.by = "sub.cluster", label=T) | FeaturePlot(vale_clean, order=T, features = genes_mesGATA2_to_neural )& scale_colour_gradientn(colours = c( "cyan", "deepskyblue", "forestgreen", "darkorange", "darkred"))  )
print(DimPlot(vale_clean, group.by = "sub.cluster", label=T) | FeaturePlot(vale_clean, order=T, features = genes_mesGATA2_to_neural2 )& scale_colour_gradientn(colours = c( "cyan", "deepskyblue", "forestgreen", "darkorange", "darkred"))  )
print(VlnPlot(vale_clean, group.by = "sub.cluster", features = genes_germcells_to_mesGATA4 ) )
print(VlnPlot(vale_clean, group.by = "sub.cluster", features = genes_mesGATA2_to_neural ) )
print(VlnPlot(vale_clean, group.by = "sub.cluster", features = genes_mesGATA2_to_neural2 ) )

# UMAP with eD vs cR
print(DimPlot(vale_clean, group.by = "comparison", sizes.highlight=0.1, reduction = "umap", label=T, cols = c("grey", "blue", "red")) )  #+scale_color_manual(labels = c( "eD", "both", "cR",), values = c("blue", "red"))    )
print(DimPlot(vale_clean, group.by = "comparison", sizes.highlight=0.1, reduction = "umap", label=F, cols = c("grey", "blue", "red")) )  #+scale_color_manual(labels = c( "eD", "both", "cR",), values = c("blue", "red"))    )


# plot venn diagram eD vs cR
listInput <- list(eD = colnames(vale_clean)[vale_clean$comparison %in% c("Emptydrops-multiome", "both")], 
                  cR = colnames(vale_clean)[vale_clean$comparison %in% c("cellRanger-arc"  , "both")]
                  )
overlaps <- euler(listInput, shape = "ellipse")
p2 <- plot(overlaps, 
           quantities = TRUE,
           labels = list(font = 4))
print(p2)

# find markers plot heatmap
eD_all.markers <- FindAllMarkers(vale_clean,  min.pct = 0.25, logfc.threshold = 0.25)
eD_all.markers50 <- eD_all.markers %>%
  group_by(cluster) %>%
  slice_max(n = 50, order_by = avg_log2FC)

eD_all.markers <- eD_all.markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)

write.table(eD_all.markers50,
            markers_tsv,
            sep = '\t', row.names = T, col.names = T, quote = F)

print(DoHeatmap(vale_clean, features = eD_all.markers$gene, size=4,
                angle = 90) + NoLegend()+ theme(axis.text.y = element_text(size = 5))   )

# calculate cR/eD venn diagram by cluster
eD_cells = colnames(vale_clean)[vale_clean$comparison %in% c("Emptydrops-multiome", "both")]
cR_cells = colnames(vale_clean)[vale_clean$comparison %in% c("cellRanger-arc"  , "both")]
venn_df <- data.frame("0" = c(as.integer(sum(vale_clean$sub.cluster=="0")),
                              as.integer(sum(vale_clean$sub.cluster=="0" & colnames(vale_clean) %in% eD_cells)),
                              as.integer(sum(vale_clean$sub.cluster=="0" & colnames(vale_clean) %in% cR_cells)),
                              sum(vale_clean$sub.cluster=="0" & colnames(vale_clean) %in% cR_cells & colnames(vale_clean) %in% eD_cells),
                              sum(vale_clean$sub.cluster=="0" & colnames(vale_clean) %in% eD_cells)/sum(vale_clean$sub.cluster=="0"),
                              sum(vale_clean$sub.cluster=="0" & colnames(vale_clean) %in% cR_cells)/sum(vale_clean$sub.cluster=="0"),
                              sum(vale_clean$sub.cluster=="0" & colnames(vale_clean) %in% cR_cells & colnames(vale_clean) %in% eD_cells)/sum(vale_clean$sub.cluster=="0")
))
for ( cl in c(1:(max(as.integer(vale_clean$sub.cluster))-1))   ){
  char_cl = as.character(cl)
  temp_df <- data.frame("new" = c(sum(vale_clean$sub.cluster==char_cl),
                                  sum(vale_clean$sub.cluster==char_cl & colnames(vale_clean) %in% eD_cells),
                                  sum(vale_clean$sub.cluster==char_cl & colnames(vale_clean) %in% cR_cells),
                                  sum(vale_clean$sub.cluster==char_cl & colnames(vale_clean) %in% cR_cells & colnames(vale_clean) %in% eD_cells),
                                  sum(vale_clean$sub.cluster==char_cl & colnames(vale_clean) %in% eD_cells)/sum(vale_clean$sub.cluster==char_cl),
                                  sum(vale_clean$sub.cluster==char_cl & colnames(vale_clean) %in% cR_cells)/sum(vale_clean$sub.cluster==char_cl),
                                  sum(vale_clean$sub.cluster==char_cl & colnames(vale_clean) %in% cR_cells & colnames(vale_clean) %in% eD_cells)/sum(vale_clean$sub.cluster==char_cl)
  ))
  venn_df <- cbind(venn_df, temp_df)
}
colnames(venn_df) <- as.character(c(0:(max(as.integer(vale_clean$sub.cluster))-1) ))
rownames(venn_df) <- c("total cells",
                       "# of eD cells",
                       "# of cR cells",
                       "# of common cells",
                       "% of eD cells",
                       "% of cR cells",
                       "% of common cells"
)
#write.table(venn_df, venns_per_cl, sep = '\t', row.names = T, col.names = T, quote = F)
eD_tpr_only <- unname(unlist(venn_df["% of eD cells",])) - unname(unlist(venn_df["% of common cells",]) )
cR_tpr_only <- unname(unlist(venn_df["% of cR cells",])) - unname(unlist(venn_df["% of common cells",]) )
common <- unname(unlist(venn_df["% of common cells",]))

Values <- matrix(c(cR_tpr_only, common, eD_tpr_only), nrow = 3, ncol = max(as.integer(vale_clean$sub.cluster)), byrow = TRUE)
barplot(Values, main = "eD vs cR by cluster", names.arg = seq(0, max(as.integer(vale_clean$sub.cluster))-1) , 
        xlab = "cluster", ylab = "fraction", col = c("grey", "pink", "salmon"))
#legend(25, 0.1, lwd=3, col=c("salmon", "pink", "grey"), lty=c(1,1,1), legend=c("EmptyDrops_multiome ", "common", "cellRanger-arc") )
legend("topleft", bg="white", lwd=3, col=c("salmon", "pink", "grey"), lty=c(1,1,1),
       legend=c("EmptyDrops_multiome ", "common", "cellRanger-arc") )


# atac clustering
srat_atac <- vale_clean
#cl4rna = colnames(vale_clean)[vale_clean$seurat_clusters=="4"]
DefaultAssay(srat_atac) <- "ATAC"
srat_atac <- RunTFIDF(srat_atac)
srat_atac <- FindTopFeatures(srat_atac, min.cutoff = 'q0')
srat_atac <- RunSVD(srat_atac)
ElbowPlot(srat_atac, n=40)
DepthCor(srat_atac, n=33)
srat_atac <- RunUMAP(object = srat_atac, reduction = 'lsi', dims = 2:30, seed.use = 6)
srat_atac <- FindNeighbors(object = srat_atac, reduction = 'lsi', dims = 2:30, random.seed = 0)
srat_atac <- FindClusters(object = srat_atac, verbose = FALSE, algorithm = 3, resolution = 1)
DefaultAssay(srat_atac) <- "SCT"
print(DimPlot(object = srat_atac, reduction = 'umap', label = TRUE)+ggtitle("clustering based on atac")) #| DimPlot(object = srat_combined, label = TRUE)+ggtitle("clustering based on rna")  
print(FeaturePlot(srat_atac, reduction = 'umap',features="nCount_RNA")& scale_colour_gradientn(colours = c( "cyan", "deepskyblue", "forestgreen", "darkorange", "darkred"),limits = c(0, 30000)) )
print(FeaturePlot(srat_atac, reduction = 'umap',order=T, features=c("DCC", "ZP3") )& scale_colour_gradientn(colours = c( "cyan", "deepskyblue", "forestgreen", "darkorange", "darkred"))  )
print(DimPlot(srat_atac,group.by = "comparison", sizes.highlight=0.1, reduction = "umap", label=T, cols = c("grey", "blue", "red")) )  #+scale_color_manual(labels = c( "eD", "both", "cR",), values = c("blue", "red"))    )
print(DimPlot(srat_atac,group.by = "comparison", sizes.highlight=0.1, reduction = "umap", label=F, cols = c("grey", "blue", "red")) )  #+scale_color_manual(labels = c( "eD", "both", "cR",), values = c("blue", "red"))    )
# in ATAC: germ cells
print(DimPlot(srat_atac, reduction = 'umap',label=T) | FeaturePlot(srat_atac, reduction = 'umap',order=T, features = genes_PGC_FGC_FGCmitotic_oogoniaSTRA8 )& scale_colour_gradientn(colours = c( "cyan", "deepskyblue", "forestgreen", "darkorange", "darkred"))  )
print(DimPlot(srat_atac, reduction = 'umap',label=T) | FeaturePlot(srat_atac, order=T, reduction = 'umap',features = genes_oogoniameiotic_preoocyte_oocyte_prespermatogonia )& scale_colour_gradientn(colours = c( "cyan", "deepskyblue", "forestgreen", "darkorange", "darkred"))  )
print(VlnPlot(srat_atac,features = genes_PGC_FGC_FGCmitotic_oogoniaSTRA8 )  )
print(VlnPlot(srat_atac,features = genes_oogoniameiotic_preoocyte_oocyte_prespermatogonia ) )
print(VlnPlot(srat_atac, features = c("DCC", "ZP3") )  )
# in ATAC: developing ovary
print(DimPlot(srat_atac, reduction = 'umap',label=T) | FeaturePlot(srat_atac, reduction = 'umap',order=T, features = genes_germcells_to_mesGATA4 )& scale_colour_gradientn(colours = c( "cyan", "deepskyblue", "forestgreen", "darkorange", "darkred"))  )
print(DimPlot(srat_atac, reduction = 'umap',label=T) | FeaturePlot(srat_atac, reduction = 'umap',order=T, features = genes_mesGATA2_to_neural )& scale_colour_gradientn(colours = c( "cyan", "deepskyblue", "forestgreen", "darkorange", "darkred"))   )
print(VlnPlot(srat_atac, features = genes_germcells_to_mesGATA4 )  )
print(VlnPlot(srat_atac, features = genes_mesGATA2_to_neural )   )
# in ATAC: qc
print(VlnPlot(srat_atac, features="percent.mt")+ stat_summary(fun.y = median, geom='point', size = 2, colour = "blue") )
print(VlnPlot(srat_atac, features="percent.ribo")+ stat_summary(fun.y = median, geom='point', size = 2, colour = "blue"))
print(VlnPlot(srat_atac, features="nFeature_RNA")+ stat_summary(fun.y = median, geom='point', size = 2, colour = "blue"))
print(VlnPlot(srat_atac, features="nCount_RNA")+ stat_summary(fun.y = median, geom='point', size = 2, colour = "blue"))
print(VlnPlot(srat_atac, features="nCount_ATAC")+ stat_summary(fun.y = median, geom='point', size = 2, colour = "blue"))
print(VlnPlot(srat_atac, features="nFeature_ATAC")+ stat_summary(fun.y = median, geom='point', size = 2, colour = "blue"))

# calculate cR/eD venn diagram by *ATAC* cluster
venn_df <- data.frame("0" = c(as.integer(sum(srat_atac$seurat_clusters=="0")),
                              as.integer(sum(srat_atac$seurat_clusters=="0" & colnames(srat_atac) %in% eD_cells)),
                              as.integer(sum(srat_atac$seurat_clusters=="0" & colnames(srat_atac) %in% cR_cells)),
                              sum(srat_atac$seurat_clusters=="0" & colnames(srat_atac) %in% cR_cells & colnames(srat_atac) %in% eD_cells),
                              sum(srat_atac$seurat_clusters=="0" & colnames(srat_atac) %in% eD_cells)/sum(srat_atac$seurat_clusters=="0"),
                              sum(srat_atac$seurat_clusters=="0" & colnames(srat_atac) %in% cR_cells)/sum(srat_atac$seurat_clusters=="0"),
                              sum(srat_atac$seurat_clusters=="0" & colnames(srat_atac) %in% cR_cells & colnames(srat_atac) %in% eD_cells)/sum(srat_atac$seurat_clusters=="0")
))
for ( cl in c(1:(max(as.integer(srat_atac$seurat_clusters))-1))   ){
  char_cl = as.character(cl)
  temp_df <- data.frame("new" = c(sum(srat_atac$seurat_clusters==char_cl),
                                  sum(srat_atac$seurat_clusters==char_cl & colnames(srat_atac) %in% eD_cells),
                                  sum(srat_atac$seurat_clusters==char_cl & colnames(srat_atac) %in% cR_cells),
                                  sum(srat_atac$seurat_clusters==char_cl & colnames(srat_atac) %in% cR_cells & colnames(srat_atac) %in% eD_cells),
                                  sum(srat_atac$seurat_clusters==char_cl & colnames(srat_atac) %in% eD_cells)/sum(srat_atac$seurat_clusters==char_cl),
                                  sum(srat_atac$seurat_clusters==char_cl & colnames(srat_atac) %in% cR_cells)/sum(srat_atac$seurat_clusters==char_cl),
                                  sum(srat_atac$seurat_clusters==char_cl & colnames(srat_atac) %in% cR_cells & colnames(srat_atac) %in% eD_cells)/sum(srat_atac$seurat_clusters==char_cl)
  ))
  venn_df <- cbind(venn_df, temp_df)
}
colnames(venn_df) <- as.character(c(0:(max(as.integer(srat_atac$seurat_clusters))-1) ))
rownames(venn_df) <- c("total cells",
                       "# of eD cells",
                       "# of cR cells",
                       "# of common cells",
                       "% of eD cells",
                       "% of cR cells",
                       "% of common cells"
)
#write.table(venn_df, venns_per_ATAC_cl, sep = '\t', row.names = T, col.names = T, quote = F)
eD_tpr_only <- unname(unlist(venn_df["% of eD cells",])) - unname(unlist(venn_df["% of common cells",]) )
cR_tpr_only <- unname(unlist(venn_df["% of cR cells",])) - unname(unlist(venn_df["% of common cells",]) )
common <- unname(unlist(venn_df["% of common cells",]))

Values <- matrix(c(cR_tpr_only, common, eD_tpr_only), nrow = 3, ncol = max(as.integer(srat_atac$seurat_clusters)), byrow = TRUE)
barplot(Values, main = "eD vs cR by cluster", names.arg = seq(0, max(as.integer(srat_atac$seurat_clusters))-1) , 
        xlab = "cluster", ylab = "fraction", col = c("grey", "pink", "salmon"))
#legend(25, 0.1, lwd=3, col=c("salmon", "pink", "grey"), lty=c(1,1,1), legend=c("EmptyDrops_multiome ", "common", "cellRanger-arc") )
legend("topleft", bg="white", lwd=3, col=c("salmon", "pink", "grey"), lty=c(1,1,1),
       legend=c("EmptyDrops_multiome ", "common", "cellRanger-arc") )

UMAP_comparison(srat_atac, vale_clean)
UMAP_comparison(vale_clean, srat_atac)

saveRDS(vale_clean, srat_vale_clean)
saveRDS(srat_atac, srat_atac_vale_clean)

# ---------------- plots FRiP-after -----------------------------------
A = vale_clean$FRiP[ colnames(vale_clean) %in% eD_cells   ]
B = vale_clean$FRiP[ colnames(vale_clean) %in% cR_cells   ]
# PLOT OUTLINE FRiP after improved plot
limits <- range(c(A,B))
breaks <- seq(limits[1], limits[2], length.out=200)
hg_eD =   hist(A , breaks=breaks, plot = FALSE)
hg_cR =   hist(B , breaks=breaks, plot = FALSE)
xrange <- range(hg_eD$mids)
yrange <- range(hg_eD$counts)
plot(0,0,type="n", xlim=xrange, ylim=yrange, xlab="", 
     ylab="", cex.axis=1.2, cex.lab=1.4, main="", font=2, cex.main=1.4)
mtext(side=1, line=2, "FRiP", col="black", font=2,cex=1.5)
mtext(side=2, line=2, "# of cells", col="black", font=2,cex=1.5)
mtext(side=3, line=0.5, sample, col="black", font=2, cex=1.5)
shift <- 0
modes = c(hg_eD$mids) #, hg_cR$mids)
plotHistogramOutline(breaks+shift, hg_cR$counts, col="grey", lwd=2)
shift <- 0
plotHistogramOutline(breaks+shift, hg_eD$counts, col="salmon", lwd=2)
legend("topright", col=c("grey","salmon"), legend=c( "cellRanger-arc", "EmptyDrops_multiome"), lwd=2, cex=1.2)
# --------------------------------------------------------

slope = -0.888259
k_means_intercept = -0.888259

# ------------------   plot combination of log counts after  -----------------
# Having a look at the distribution of retained cells in more detail.
e.keep <- eD.out_multi$Row.names %in% eD_cells
c.keep <- eD.out_multi$Row.names %in% cR_cells
union_cells = union(eD_cells, cR_cells)
log_comb <- log10(eD.out_multi$Total_RNA) - slope * log10(eD.out_multi$Total_chromatin)
limits <- range(log_comb[eD.out_multi$Row.names %in% union_cells])
breaks <- seq(limits[1], limits[2], length.out=21)
modes <- list("EmptyDrops"=e.keep, "CellRanger"=c.keep)

collected.x <- collected.y <- list()
for (mode in names(modes)) {
  keep <- modes[[mode]]
  d <- hist(log_comb[keep], breaks=breaks, plot=FALSE)
  collected.x[[mode]] <-d$mid
  collected.y[[mode]] <-d$count
}
xrange <- range(sapply(collected.x, range))
yrange <- range(sapply(collected.y, range))

plot(0,0,type="n", xlim=xrange, ylim=yrange, xlab=expression(Log[10]~"(ATAC count)"-"slope * "~Log[10]~"(RNA count)"), 
     ylab="Number of cells", cex.axis=1.2, cex.lab=1.4, main=sample, cex.main=1.4)

shift <- 0
colors <- c("EmptyDrops"="salmon", "CellRanger"="grey")
for (mode in names(modes)) { 
  plotHistogramOutline(breaks+shift, collected.y[[mode]], col=colors[mode], lwd=2)
  shift <- shift + 0.01
}

legend("topright", col=colors[names(modes)], legend=names(modes), lwd=2, cex=1.2)
abline(v=k_means_intercept)

hist(log_comb, breaks=200, xlab=expression(Log[10]~"(ATAC count)"-"slope * "~Log[10]~"(RNA count)"))
rect(0, 0, k_means_intercept, par("usr")[4], col = "pink", border = NA)
par(new = TRUE)
rect(k_means_intercept, 0 ,9, par("usr")[4], col = "lightblue", border = NA)
par(new = TRUE)
hist(log_comb, ylim= c(0,5000), breaks=500)
abline(v=k_means_intercept, col="red")
# ---------------------------------------------------------------

dev.off()




}
