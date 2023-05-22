setwd("/mnt/beegfs6/home3/ahringer/em613/analysis/multiomics/emptryDrops_multiome2023_eDv3")
library(Seurat)
library(Signac)
library(DropletUtils)
library(ggplot2)
library(ggridges)
library(mixtools)
library(UpSetR)
library(eulerr)
library(DelayedMatrixStats)
library(eDv3)
library(dplyr)
#library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
source("real/vale_annotate.R")
source("real/clustering_comparison.R")
source("simulations/fcn_for_sim.R")
current_date= paste(unlist(strsplit(as.character(Sys.Date()), "-")), collapse="")
opath <- paste0("data/output/realdata/", current_date)
dir.create(file.path(opath), recursive=TRUE)


ALLFILES <- c(#"pbmc_gran_sorted/raw_feature_bc_matrix.h5",
              # "gastr_d3.5_A/gastr_d3_5_multiome_A_L001_raw_feature_bc_matrix.h5",
              # "gastr_d3.5_B/gastr_d3_5_multiome_B_L001_raw_feature_bc_matrix.h5",
              # "gastr_d4.5_A/gastr_d4_5_multiome_A_L001_raw_feature_bc_matrix.h5",
              # "gastr_d4.5_B/gastr_d4_5_multiome_B_L001_raw_feature_bc_matrix.h5",
              # "gastr_d5_A/gastr_d5_multiome_A_L001_raw_feature_bc_matrix.h5",
              # "gastr_d5_B/gastr_d5_multiome_B_L001_raw_feature_bc_matrix.h5",
              # 
              # "gastr_d4_A/raw_feature_bc_matrix.h5",
              # "gastr_d4_B/gastr_d4_multiome_B_L001_raw_feature_bc_matrix.h5"
               "valentina_8176/FCA_GND10288176_raw_feature_bc_matrix.h5",
               "valentina_8177/FCA_GND10288177_raw_feature_bc_matrix.h5"
              # "valentina_8178/FCA_GND10288178_raw_feature_bc_matrix.h5",
              # "valentina_8179/FCA_GND10288179_raw_feature_bc_matrix.h5",
              # "valentina_8180/FCA_GND10288180_raw_feature_bc_matrix.h5"
              
            )

ALL_BARCODES <- c(#"pbmc_gran_sorted/barcodes.tsv.gz",
                  # "gastr_d3.5_A/gastr_d3_5_multiome_A_L001_barcodes.tsv.gz",
                  # "gastr_d3.5_B/gastr_d3_5_multiome_B_L001_barcodes.tsv.gz",
                  # "gastr_d4.5_A/gastr_d4_5_multiome_A_L001_barcodes.tsv.gz",
                  # "gastr_d4.5_B/gastr_d4_5_multiome_B_L001_barcodes.tsv.gz",
                  # "gastr_d5_A/gastr_d5_multiome_A_L001_barcodes.tsv.gz",
                  # "gastr_d5_B/gastr_d5_multiome_B_L001_barcodes.tsv.gz",
                  # 
                  # "gastr_d4_A/gastr_d4_multiome_A_L001_barcodes.tsv.gz",
                  # "gastr_d4_B/gastr_d4_multiome_B_L001_barcodes.tsv.gz"
                  "valentina_8176/barcodes.tsv.gz",
                  "valentina_8177/barcodes.tsv.gz"
                  # "valentina_8178/barcodes.tsv.gz",
                  # "valentina_8179/barcodes.tsv.gz",
                  # "valentina_8180/barcodes.tsv.gz"
)


for (i in seq_along(ALLFILES) ) { 

fname = ALLFILES[i]
sce <- Read10X_h5(file.path("data/input", fname))

stub <- sub("/.*", "", fname, "_qc")
cR_barcodes_file <-  gzfile(file.path("data/input", ALL_BARCODES[i] ))

if (substring(stub, 1,4)=="vale"){
metadata10x <- read.csv(
  file = paste0("data/input/", stub, "/FCA_GND1028",  substring(stub, 11,14), "_per_barcode_metrics.csv"),
  header = TRUE,
  row.names = 1
)
}

dir.create(file.path(opath, stub), recursive=TRUE)


#ofile <- file.path(opath, paste0(stub, "_eD-multiome-res.tsv"))
ffile <- file.path(opath, stub, paste0(stub, "_eD_multiome.pdf"))  
ffile2 <- file.path(opath, stub, paste0(stub, "_compare_to_cR.pdf"))
ffile_qc <- file.path(opath, stub, paste0(stub, "_ribo_mito_qc.pdf"))
f_cts_log <- file.path(opath, stub, paste0(stub, "_hist_counts_log_space.pdf"))  
heatfile <- file.path(opath, stub, paste0(stub, "_cluster_heatmap.pdf"))  
eD_multi_tsv <- file.path(opath, stub, paste0(stub, "_eD_multiome.tsv")) 
markers_tsv <- file.path(opath, stub, paste0(stub, "_cluster_markers.tsv")) 
srat_file <- file.path(opath, stub, paste0("srat_", stub, ".rds"))  
srat_atac_file <- file.path(opath, stub, paste0("srat_atac_", stub, ".rds"))  
corr_csv <- file.path(opath, stub, paste0("corr_", stub, ".csv"))  
corr_inv_csv <- file.path(opath, stub, paste0("corr_inv_", stub, ".csv"))  
corr_csv_fine <- file.path(opath, stub, paste0("corr_fine_", stub, ".csv"))  
corr_inv_csv_fine <- file.path(opath, stub, paste0("corr_inv_fine_", stub, ".csv"))  
venns_per_cl <- file.path(opath, stub, paste0("venns_per_cl_", stub, ".csv"))  
venns_per_cl_pdf <- file.path(opath, stub, paste0("venns_per_cl_", stub, ".pdf"))  
cl_res1 <- file.path(opath, stub, paste0("clustering_res1_", stub, ".pdf"))  
ATAC_clust <- file.path(opath, stub, paste0("ATAC_clustering_", stub, ".pdf"))  
venns_per_ATAC_cl <- file.path(opath, stub, paste0("venns_per_ATAC_cl_", stub, ".csv"))  
srat_cR_file <- file.path(opath, stub, paste0("srat_cR_", stub, ".pdf"))  
eD_metadata <- file.path(opath, stub, paste0("eD_metadata_", stub, ".tsv"))  

count_matrix_rna <- sce[["Gene Expression"]]
count_matrix_atac <- sce[["Peaks"]]

pdf(ffile)


lower_rna = NULL
barhop_rna = NULL
lower_atac = NULL
barhop_atac = NULL

start_time <- Sys.time()
set.seed(0)
eD.out_multi <- emptydrops_multiome(count_matrix_rna, lower_rna, barhop_rna, count_matrix_atac, lower_atac, barhop_atac )
print("the number of cells detected is: ")
print(sum(eD.out_multi$FDR_multi<0.001 & ! is.na(eD.out_multi$FDR_multi)))
end_time <- Sys.time()
print(end_time - start_time)


dev.off()


write.table(eD.out_multi,
            paste0(eD_multi_tsv),
            sep = '\t', row.names = T, col.names = T, quote = F)

write.table(eD.out_multi@metadata,
            eD_metadata,
            sep = '\t', row.names = T, col.names = T, quote = F)





cR_cells <- readLines(cR_barcodes_file)

if (substring(stub, 1,4)=="vale"){
  srat <- CreateSeuratObject(counts = count_matrix_rna,
                           meta.data = metadata10x)
  srat$FRiP <- srat$atac_peak_region_fragments / srat$atac_fragments * 100
  
  srat[["ATAC"]] <- CreateChromatinAssay(
    counts = count_matrix_atac,
    min.features = -1,
    sep = c(":", "-"),
    genome = "hg38",
    fragments = paste0("data/input/", stub, "/FCA_GND1028",  substring(stub, 11,14), "_atac_fragments.tsv.gz") ,
    #  annotation = genomeAnnotation$chromSizes,
    validate.fragments = FALSE
  )
  
  srat$FDR_multi <- eD.out_multi[colnames(srat),]$FDR_multi
  srat$k_means <- eD.out_multi[colnames(srat),]$k_means
  srat$FDR <- eD.out_multi[colnames(srat),]$FDR
  
  
  
}else{
  srat <- CreateSeuratObject(counts = count_matrix_rna)
  
  srat[["ATAC"]] <- CreateChromatinAssay(
    counts = count_matrix_atac,
    min.features = -1,
    sep = c(":", "-"),
    #fragments = paste0("data/input/", stub, "/FCA_GND1028",  substring(stub, 11,14), "_atac_fragments.tsv.gz") ,
    #  annotation = genomeAnnotation$chromSizes,
    validate.fragments = FALSE
  )
}



if (substring(stub, 1,4)=="vale"){
  
  # compute nucleosome signal score per cell
  DefaultAssay(srat) <- "ATAC"
  srat <- NucleosomeSignal(object = srat)

  # extract gene annotations from EnsDb
  annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
  
  # change to UCSC style
  seqlevelsStyle(annotations) <- 'UCSC'
  
  # add the gene information to the object
  Annotation(srat) <- annotations
  
  # compute TSS enrichment score per cell
  srat <- TSSEnrichment(object = srat, fast = FALSE)
}
DefaultAssay(srat) <- "RNA"

pdf(ffile2)

print(cR_barcodes_file)
cR_cells <- readLines(cR_barcodes_file)
eD_cells <- eD.out_multi$Row.names[ eD.out_multi$FDR_multi<0.001 & ! is.na(eD.out_multi$FDR_multi) ]
eD_rna_cells <- eD.out_multi$Row.names[ eD.out_multi$FDR_RNA<0.001 & ! is.na(eD.out_multi$FDR_RNA) ]
listInput <- list(eD = eD_cells, 
                  cR = cR_cells,
                  eDrna =eD_rna_cells
                  )
p1 <- upset(fromList(listInput), nsets = 6,, order.by = "freq")
overlaps <- euler(listInput, shape = "ellipse")
p2 <- plot(overlaps, 
           quantities = TRUE,
           labels = list(font = 4))
print(p1)
print(p2)

listInput <- list(eD = eD_cells, 
                  cR = cR_cells
                  )
p1 <- upset(fromList(listInput), nsets = 6,, order.by = "freq")
overlaps <- euler(listInput, shape = "ellipse")
p2 <- plot(overlaps, 
           quantities = TRUE,
           labels = list(font = 4))
print(p1)
print(p2)

# make umap
below_k_means <- rownames(eD.out_multi)[!as.logical(eD.out_multi$k_means)]
bkmeans_eD <- intersect(eD_cells, below_k_means)
akmeans_eD <- setdiff(eD_cells, below_k_means)
union_cells <- union(eD_cells, cR_cells)
eD_minus_cR <- setdiff(eD_cells, cR_cells)
cR_minus_eD <- setdiff(cR_cells, eD_cells )
intersection <- intersect(eD_cells, cR_cells)
srat_subset <- subset(srat, cells=union_cells )

# calculate QC metrics: mito and ribo contents and FRiP
if (substring(stub, 1,4)=="vale"){
  
  FRiP <- srat$atac_peak_region_fragments / srat$atac_fragments * 100
  hist(FRiP, breaks=500)
  #abline(v=3)
  abline(v=0.5)
  hist(FRiP, breaks=500, ylim = c(0,100))
  hist(FRiP[srat$is_cell==1] , breaks=500, ylim = c(0,100))
  #hist(named_FRiP[eD_cells] , breaks=500, ylim = c(0,100))
  
    
  C <- srat_subset@assays[["RNA"]]
  rb.genes <- rownames(C)[grep("^RP[SL]",rownames(C))]
  percent.ribo <- colSums(C[rb.genes,])/Matrix::colSums(C)*100
  srat_subset <- AddMetaData(srat_subset, percent.ribo, col.name = "percent.ribo")
  ribo_lim <- 100.1  #median(srat_subset[["percent.ribo"]][,1]) + 3* mad(srat_subset[["percent.ribo"]][,1])
  
  srat_subset[["percent.mt"]] <- PercentageFeatureSet(srat_subset, pattern = "^MT-")
  mito_lim <- median(srat_subset[["percent.mt"]][,1]) + 3* mad(srat_subset[["percent.mt"]][,1])

  }else{
  
  C <- srat_subset@assays[["RNA"]]
  rb.genes <- rownames(C)[grep("^Rp[sl]",rownames(C))]
  percent.ribo <- colSums(C[rb.genes,])/Matrix::colSums(C)*100
  srat_subset <- AddMetaData(srat_subset, percent.ribo, col.name = "percent.ribo")
  ribo_lim <- median(srat_subset[["percent.ribo"]][,1]) + 3* mad(srat_subset[["percent.ribo"]][,1])
  
  srat_subset[["percent.mt"]] <- PercentageFeatureSet(srat_subset, pattern = "^mt-")
  mito_lim <- median(srat_subset[["percent.mt"]][,1]) + 3* mad(srat_subset[["percent.mt"]][,1])
    
}


hist(srat_subset[["percent.ribo"]][,1] , breaks=100)
abline(v=ribo_lim)
hist(srat_subset[["percent.mt"]][,1] , breaks=100)
abline(v=mito_lim)



if (substring(stub, 1,4)=="vale"){
  
  hist(srat_subset[["TSS.enrichment"]][,1] , breaks=100)
  
  print(VlnPlot(
    object = srat_subset,
    features = c(
      'TSS.enrichment'),
    pt.size = 0.1
  )+NoLegend() )
  
  hist(srat$FRiP, breaks=500)
  abline(v=3)
  abline(v=0.5)
  
  hist(srat$FRiP, breaks=500, ylim=c(0,50000))
  abline(v=3)
  abline(v=0.5)
  
  hist(srat_subset$FRiP, breaks=500)
  abline(v=3)
  abline(v=0.5)
}

plot.new()
retained = sum(srat_subset[["percent.ribo"]][,1] < ribo_lim & srat_subset[["percent.mt"]][,1] < mito_lim )
retained_in_cR = colnames(srat_subset)[srat_subset[["percent.ribo"]][,1] < ribo_lim & srat_subset[["percent.mt"]][,1] < mito_lim & colnames(srat_subset) %in% cR_cells]
retained_in_eD = colnames(srat_subset)[srat_subset[["percent.ribo"]][,1] < ribo_lim & srat_subset[["percent.mt"]][,1] < mito_lim & colnames(srat_subset) %in% eD_cells]
rejected = sum(srat_subset[["percent.ribo"]][,1] > ribo_lim | srat_subset[["percent.mt"]][,1] > mito_lim )
text(x=0.2, y=.1, paste0("retained=", retained))
text(x=0.2, y=0.2, paste0("rejected=", rejected))
text(x=0.2, y=0.3, paste0("retained_in_cR=", length(retained_in_cR) ) )
text(x=0.2, y=0.4, paste0("retained_in_eD=", length(retained_in_eD)) )
text(x=0.2, y=0.5, paste0("slope=", eD.out_multi@metadata[["k_means_slope"]]))
text(x=0.2, y=0.6, paste0("intercept=", eD.out_multi@metadata[["k_means_intercept"]] ))

listInput <- list(retained_eD = retained_in_eD, 
                  retained_cR = retained_in_cR
)
p1 <- upset(fromList(listInput), nsets = 6,, order.by = "freq")
overlaps <- euler(listInput, shape = "ellipse")
p2 <- plot(overlaps, 
           quantities = TRUE,
           labels = list(font = 4))
print(p1)
print(p2)



#srat_subset <- subset(srat, cells=union_cells )
srat_subset <- subset(x = srat_subset, subset = percent.mt < mito_lim & percent.ribo < ribo_lim)
srat_subset <- SCTransform(srat_subset)
srat_subset <- RunPCA(srat_subset, seed.use=42, features = VariableFeatures(object = srat_subset))
print(ElbowPlot(srat_subset, ndims = 50)     )
srat_subset <- FindNeighbors(srat_subset, dims = 1:50)
srat_subset <- FindClusters(srat_subset, resolution = 2, random.seed = 0)
srat_subset <- RunUMAP(srat_subset, dims = 1:50, seed.use = 42)

srat_subset$comparison <- 1
srat_subset$comparison[ colnames(srat_subset) %in% eD_minus_cR] <- "Emptydrops-multiome"
srat_subset$comparison[ colnames(srat_subset) %in% intersection] <- "both"
srat_subset$comparison[ colnames(srat_subset) %in% cR_minus_eD] <- "cellRanger-arc"
print(DimPlot(srat_subset, group.by = "comparison", sizes.highlight=0.1, reduction = "umap", label=T, cols = c("grey", "blue", "red")) )  #+scale_color_manual(labels = c( "eD", "both", "cR",), values = c("blue", "red"))    )

print(DimPlot(srat_subset, reduction = "umap", label=T)     )
print(DimPlot(srat_subset, cells.highlight = eD_minus_cR, sizes.highlight=0.1, reduction = "umap", label=T)+scale_color_manual(labels = c("other", "eD_minus_cR"), values = c("blue", "red"))    )
print(DimPlot(srat_subset, cells.highlight = intersection, sizes.highlight=0.1, reduction = "umap", label=T)+scale_color_manual(labels = c("other", "intersection"), values = c("blue", "red"))  )
print(DimPlot(srat_subset, cells.highlight = eD_minus_cR, sizes.highlight=0.1, reduction = "umap", label=T)+scale_color_manual(labels = c("other", "eD_minus_cR"), values = c("grey", "red"))    )
print(DimPlot(srat_subset, cells.highlight = cR_minus_eD, sizes.highlight=0.1, reduction = "umap", label=T)+scale_color_manual(labels = c("other", "cR_minus_eD"), values = c("grey", "red"))    )
print(DimPlot(srat_subset, cells.highlight = bkmeans_eD, sizes.highlight=0.1, reduction = "umap", label=T)+scale_color_manual(labels = c("other", "bkmeans_eD"), values = c("grey", "red"))      )
print(DimPlot(srat_subset, cells.highlight = bkmeans_eD, sizes.highlight=0.1, reduction = "umap", label=T)+scale_color_manual(labels = c("other", "bkmeans_eD"), values = c("blue", "red"))      )
print(DimPlot(srat_subset, cells.highlight = akmeans_eD, sizes.highlight=0.1, reduction = "umap", label=T)+scale_color_manual(labels = c("other", "akmeans_eD"), values = c("grey", "red"))      )
print(FeaturePlot(srat_subset, order=T, features = "nCount_RNA", ) & scale_colour_gradientn(colours = c( "cyan", "deepskyblue", "forestgreen", "darkorange", "darkred"), limits = c(1, 14000))   )
print(FeaturePlot(srat_subset, order=T, features = "percent.ribo", ) & scale_colour_gradientn(colours = c( "cyan", "deepskyblue", "forestgreen", "darkorange", "darkred")   ) )
print(FeaturePlot(srat_subset, order=T, features = "percent.mt", ) & scale_colour_gradientn(colours = c( "cyan", "deepskyblue", "forestgreen", "darkorange", "darkred"))   )  
print(FeaturePlot(srat_subset, order=T, features = "nCount_ATAC", ) & scale_colour_gradientn(colours = c( "cyan", "deepskyblue", "forestgreen", "darkorange", "darkred"), limits = c(1, 4000)) )
# FeaturePlot(srat_subset, order=T, features = c("lsl-1", "vet-6", "pes-10") ) & scale_colour_gradientn(colours = c("grey90", "cyan", "deepskyblue", "forestgreen", "darkorange", "darkred"))
# FeaturePlot(srat_subset, order=T, features = c("med-1", "med-2", "end-1", "end-3", "elt-2", "elt-7") ) & scale_colour_gradientn(colours = c("grey90", "cyan", "deepskyblue", "forestgreen", "darkorange", "darkred"))
# FeaturePlot(srat_subset, order=T, features = c("pha-4", "tbx-35", "tbx-38") ) & scale_colour_gradientn(colours = c("grey90", "cyan", "deepskyblue", "forestgreen", "darkorange", "darkred"))
print(VlnPlot(srat_subset, features = "nCount_RNA",  y.max = 13000))
print(VlnPlot(srat_subset, features = "nFeature_RNA",  y.max = 10000))
print(VlnPlot(srat_subset, features = "percent.ribo"))
print(VlnPlot(srat_subset, features = "percent.mt") + stat_summary(fun.y = median, geom='point', size = 10, colour = "blue")    )
print(VlnPlot(srat_subset, features = "nCount_ATAC",  y.max = 1000))
if (substring(stub, 1,4)=="vale"){
  FRiP_sub <- srat_subset$atac_peak_region_fragments / srat_subset$atac_fragments * 100
  max_frip_of_excluded = max(srat_subset$FRiP[srat_subset$excluded_reason==2]  )
  min_frip_of_cells = min(srat_subset$FRiP[srat_subset$is_cell==1]  )
  
  hist(FRiP_sub, breaks=500)
  abline(v=3)
  print(FeaturePlot(srat_subset, order=T, features = c("DCC", "ZP3") ) & scale_colour_gradientn(colours = c("grey90", "cyan", "deepskyblue", "forestgreen", "darkorange", "darkred")) )+ ggtitle('oocyte DCC+' )
  print(VlnPlot(srat_subset, features = c("DCC", "ZP3") ) ) 
  print(FeaturePlot(srat_subset, order=T, features = c("GDF9", "FIGLA", "NOBOX", "ZP3", "LHX8", "OOSP2") ) & scale_colour_gradientn(colours = c("grey90", "cyan", "deepskyblue", "forestgreen", "darkorange", "darkred")) )+ ggtitle('oocyte' )
  print(FeaturePlot(srat_subset, order=T, features = c("HBA1", "HBA2", "HBG1", "HBG2") ) & scale_colour_gradientn(colours = c("grey90", "cyan", "deepskyblue", "forestgreen", "darkorange", "darkred")) )+ ggtitle(' redblood ' )
  print(FeaturePlot(srat_subset, order=T, features = c("IFITM1", "NANOG", "NANOS3", "PU5F1") ) & scale_colour_gradientn(colours = c("grey90", "cyan", "deepskyblue", "forestgreen", "darkorange", "darkred")) )+ ggtitle('PGC vs FGC' )
  print(FeaturePlot(srat_subset, order=T, features = c("MEIOC", "SYCP1") ) & scale_colour_gradientn(colours = c("grey90", "cyan", "deepskyblue", "forestgreen", "darkorange", "darkred") ) )+ ggtitle('oogonia meiotic' )
  print(FeaturePlot(srat_subset, order=T, features = c("STRA8", "ZGLP1") ) & scale_colour_gradientn(colours = c("grey90", "cyan", "deepskyblue", "forestgreen", "darkorange", "darkred") ) )+ ggtitle('oogonia STRA8' )
  print(FeaturePlot(srat_subset, order=T, features = c("SPOCD1", "MORC1") ) & scale_colour_gradientn(colours = c("grey90", "cyan", "deepskyblue", "forestgreen", "darkorange", "darkred") ) )+ ggtitle('pre-spermatogonia' )

  print(VlnPlot(srat_subset, features = c("GDF9", "FIGLA", "NOBOX", "ZP3", "LHX8", "OOSP2")))
  print(VlnPlot(srat_subset, features = c("HBA1", "HBA2", "HBG1", "HBG2")))
  print(VlnPlot(srat_subset, features = c("IFITM1", "NANOG", "NANOS3", "PU5F1")   ))
  print(VlnPlot(srat_subset, features = c("MEIOC", "SYCP1")   ))
  print(VlnPlot(srat_subset, features = c("STRA8", "ZGLP1")   ))
  print(VlnPlot(srat_subset, features = c("SPOCD1", "MORC1")   ))
  print(VlnPlot(srat_subset, features = c("GATA4")   ))
  print(VlnPlot(srat_subset, features = c("PDGFRA")   ))
  print(VlnPlot(srat_subset, features = c("LGR5", "TSPAN8")   ))
  print(VlnPlot(srat_subset, features = c("PAX8", 'SOX9')   ))
  # DE 4 vs 16 in 8176
  print(VlnPlot(srat_subset, features = c("NXPH1", 'LINC02715', "ATP10B", "CNTN6", "GRM5", "DSCAM", "DCC", "SOX5", "NTN1")   ))
  # DE 14 vs 3 in 8176
  #KHDRBS2 TTC8
  
  print(VlnPlot(srat_subset, features = c('LRRN4', 'UPK3B', 'GATA4' ,  
                                          "LHX9" ,"FOXL2", "LGR5" ,"TSPAN8", 
                                          "OSR1" ,"PLAC1" ,"CYP19A1", "CALB2", 
                                          'IRX3', 'LHX2' ,'CYP26B1', 'BMP2', 'NOTCH3')   ))
  
  # germ cells
  genes_PGC_FGC_FGCmitotic_oogoniaSTRA8 = c("DAZL", "IFITM1", "NANOG", "NANOS3", "POU5F1", "DDX4", "MAEL", "ZGLP1", "STRA8")
  genes_oogoniameiotic_preoocyte_oocyte_prespermatogonia = c( "MEIOC", "SYCP1", "FIGLA", "LHX8", "NOBOX", "ZP3", "GDF9", "SPOCD1", "MORC1")
  print(DimPlot(srat_subset, label=T) | FeaturePlot(srat_subset, order=T, features = genes_PGC_FGC_FGCmitotic_oogoniaSTRA8 )& scale_colour_gradientn(colours = c( "cyan", "deepskyblue", "forestgreen", "darkorange", "darkred")) )
  print(DimPlot(srat_subset, label=T) | FeaturePlot(srat_subset, order=T, features = genes_oogoniameiotic_preoocyte_oocyte_prespermatogonia )& scale_colour_gradientn(colours = c( "cyan", "deepskyblue", "forestgreen", "darkorange", "darkred")))
  print(VlnPlot(srat_subset, features = genes_PGC_FGC_FGCmitotic_oogoniaSTRA8 ) )
  print(VlnPlot(srat_subset, features = genes_oogoniameiotic_preoocyte_oocyte_prespermatogonia ))
  print(VlnPlot(srat_subset, features = c("DCC", "ZP3") ))
  
  # developing ovary
  genes_germcells_to_mesGATA4 = c( "DAZL", "UPK3B", "GATA4", "LHX9", "NR5A1", "WNT6", "IRX3", "FOXL2", "ARX")
  genes_mesGATA2_to_neural = c( "TCF21", "PDGFRA", "DCV", "GATA2", "NR2F1", "PDGFRB")
  genes_mesGATA2_to_neural2 = c( "MYH11", "PTPRC", "CDH5", "PAX8", "EPCAM", "HBA1", "ASCL1")
  print(DimPlot(srat_subset, label=T) | FeaturePlot(srat_subset, order=T, features = genes_germcells_to_mesGATA4 )& scale_colour_gradientn(colours = c( "cyan", "deepskyblue", "forestgreen", "darkorange", "darkred")) )
  print(DimPlot(srat_subset, label=T) | FeaturePlot(srat_subset, order=T, features = genes_mesGATA2_to_neural )& scale_colour_gradientn(colours = c( "cyan", "deepskyblue", "forestgreen", "darkorange", "darkred"))  )
  print(DimPlot(srat_subset, label=T) | FeaturePlot(srat_subset, order=T, features = genes_mesGATA2_to_neural2 )& scale_colour_gradientn(colours = c( "cyan", "deepskyblue", "forestgreen", "darkorange", "darkred"))  )
  print(VlnPlot(srat_subset, features = genes_germcells_to_mesGATA4 ) )
  print(VlnPlot(srat_subset, features = genes_mesGATA2_to_neural ) )
  print(VlnPlot(srat_subset, features = genes_mesGATA2_to_neural2 ) )
  
  
  print(RidgePlot(srat_subset, features = "FRiP" )+
          NoLegend() +
          scale_discrete_manual(aesthetics = "point_shape", values = c(21, 22, 23))+
          geom_vline(xintercept = max_frip_of_excluded)) #, ncol = 2))
  
  df_FRiP <- data.frame("frip"=srat_subset$FRiP, "cluster"=srat_subset$seurat_clusters)
  print(ggplot(df_FRiP, aes(x = frip, y = cluster, height = stat(density))) + 
    geom_density_ridges(stat = "binline", bins = 80, scale = 0.95, draw_baseline = FALSE)+
    theme_ridges(grid = FALSE, center_axis_labels = TRUE)+
    geom_vline(xintercept = max_frip_of_excluded))
  
   
  
  }else{
  print(FeaturePlot(srat_subset, order=T, features = c(     "Dppa5a", "Utf1", "Zfp42", "Pou5f1", "Brachyury" , "Wnt3"     ) ) & scale_colour_gradientn(colours = c("grey90", "cyan", "deepskyblue", "forestgreen", "darkorange", "darkred")) )
  print(FeaturePlot(srat_subset, order=T, features = c(     "Thy1", "Hes3", "Nrp2", "Epha5", "Gbx2", "Sfrp1", "Ncam1", "Pou5f1", "Hoxc8", "Hoxb9", "Hoxaas3", "Pim2", "T",  "Wnt3"     ) ) & scale_colour_gradientn(colours = c("grey90", "cyan", "deepskyblue", "forestgreen", "darkorange", "darkred")) )
  

}
  
dev.off()

saveRDS(srat_subset, srat_file)

# calculate cR/eD venn diagram by cluster
venn_df <- data.frame("0" = c(as.integer(sum(srat_subset$seurat_clusters=="0")),
                              as.integer(sum(srat_subset$seurat_clusters=="0" & colnames(srat_subset) %in% eD_cells)),
                              as.integer(sum(srat_subset$seurat_clusters=="0" & colnames(srat_subset) %in% cR_cells)),
                              sum(srat_subset$seurat_clusters=="0" & colnames(srat_subset) %in% cR_cells & colnames(srat_subset) %in% eD_cells),
                              sum(srat_subset$seurat_clusters=="0" & colnames(srat_subset) %in% eD_cells)/sum(srat_subset$seurat_clusters=="0"),
                              sum(srat_subset$seurat_clusters=="0" & colnames(srat_subset) %in% cR_cells)/sum(srat_subset$seurat_clusters=="0"),
                              sum(srat_subset$seurat_clusters=="0" & colnames(srat_subset) %in% cR_cells & colnames(srat_subset) %in% eD_cells)/sum(srat_subset$seurat_clusters=="0")
                                  ))
for ( cl in c(1:(max(as.integer(srat_subset$seurat_clusters))-1))   ){
  char_cl = as.character(cl)
  temp_df <- data.frame("new" = c(sum(srat_subset$seurat_clusters==char_cl),
                                           sum(srat_subset$seurat_clusters==char_cl & colnames(srat_subset) %in% eD_cells),
                                           sum(srat_subset$seurat_clusters==char_cl & colnames(srat_subset) %in% cR_cells),
                                           sum(srat_subset$seurat_clusters==char_cl & colnames(srat_subset) %in% cR_cells & colnames(srat_subset) %in% eD_cells),
                                           sum(srat_subset$seurat_clusters==char_cl & colnames(srat_subset) %in% eD_cells)/sum(srat_subset$seurat_clusters==char_cl),
                                           sum(srat_subset$seurat_clusters==char_cl & colnames(srat_subset) %in% cR_cells)/sum(srat_subset$seurat_clusters==char_cl),
                                           sum(srat_subset$seurat_clusters==char_cl & colnames(srat_subset) %in% cR_cells & colnames(srat_subset) %in% eD_cells)/sum(srat_subset$seurat_clusters==char_cl)
  ))
  venn_df <- cbind(venn_df, temp_df)
}
colnames(venn_df) <- as.character(c(0:(max(as.integer(srat_subset$seurat_clusters))-1) ))
rownames(venn_df) <- c("total cells",
                       "# of eD cells",
                       "# of cR cells",
                       "# of common cells",
                       "% of eD cells",
                       "% of cR cells",
                       "% of common cells"
)
write.table(venn_df, venns_per_cl, sep = '\t', row.names = T, col.names = T, quote = F)
eD_tpr_only <- unname(unlist(venn_df["% of eD cells",])) - unname(unlist(venn_df["% of common cells",]) )
cR_tpr_only <- unname(unlist(venn_df["% of cR cells",])) - unname(unlist(venn_df["% of common cells",]) )
common <- unname(unlist(venn_df["% of common cells",]))

pdf(venns_per_cl_pdf)
Values <- matrix(c(cR_tpr_only, common, eD_tpr_only), nrow = 3, ncol = max(as.integer(srat_subset$seurat_clusters)), byrow = TRUE)
barplot(Values, main = "eD vs cR by cluster", names.arg = seq(0, max(as.integer(srat_subset$seurat_clusters))-1) , 
        xlab = "cluster", ylab = "fraction", col = c("grey", "pink", "salmon"))
#legend(25, 0.1, lwd=3, col=c("salmon", "pink", "grey"), lty=c(1,1,1), legend=c("EmptyDrops_multiome ", "common", "cellRanger-arc") )
legend("topleft", bg="white", lwd=3, col=c("salmon", "pink", "grey"), lty=c(1,1,1),
       legend=c("EmptyDrops_multiome ", "common", "cellRanger-arc") )
dev.off()


# what would happens with just cR cells
pdf(srat_cR_file)
srat_cR <- subset(srat_subset, cells = cR_cells)
srat_cR <- SCTransform(srat_cR)
srat_cR <- RunPCA(srat_cR, seed.use=42, features = VariableFeatures(object = srat_cR))
print(ElbowPlot(srat_cR, ndims = 50)     )
srat_cR <- FindNeighbors(srat_cR, dims = 1:50)
srat_cR <- FindClusters(srat_cR, resolution = 2, random.seed = 0)
srat_cR <- RunUMAP(srat_cR, dims = 1:50, seed.use = 42)
print(DimPlot(srat_cR, reduction = "umap", label=T)     )
print(FeaturePlot(srat_cR, order=T, features = c("DCC", "ZP3") ) & scale_colour_gradientn(colours = c("grey90", "cyan", "deepskyblue", "forestgreen", "darkorange", "darkred") ) )+ ggtitle('UMAP with cR cells' )
dev.off()



eD_all.markers <- FindAllMarkers(srat_subset,  min.pct = 0.25, logfc.threshold = 0.25)
eD_all.markers30 <- eD_all.markers %>%
  group_by(cluster) %>%
  slice_max(n = 30, order_by = avg_log2FC)

eD_all.markers <- eD_all.markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)

write.table(eD_all.markers30,
            markers_tsv,
            sep = '\t', row.names = T, col.names = T, quote = F)

pdf(heatfile)
print(DoHeatmap(srat_subset, features = eD_all.markers$gene, size=4,
          angle = 90) + NoLegend()+ theme(axis.text.y = element_text(size = 5))   )
dev.off() 



# Having a look at the distribution of retained cells in more detail.
e.keep <- eD.out_multi$Row.names %in% eD_cells
c.keep <- eD.out_multi$Row.names %in% cR_cells
log_comb <- log10(eD.out_multi$Total_RNA) - eD.out_multi@metadata[["k_means_slope"]] * log10(eD.out_multi$Total_chromatin)
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

pdf(f_cts_log)
plot(0,0,type="n", xlim=xrange, ylim=yrange, xlab=expression(Log[10]~"(ATAC count)"-"slope * "~Log[10]~"(RNA count)"), 
     ylab="Number of cells", cex.axis=1.2, cex.lab=1.4, main=stub, cex.main=1.4)

shift <- 0
colors <- c("EmptyDrops"="salmon", "CellRanger"="grey")
for (mode in names(modes)) { 
  plotHistogramOutline(breaks+shift, collected.y[[mode]], col=colors[mode], lwd=2)
  shift <- shift + 0.01
}

legend("topright", col=colors[names(modes)], legend=names(modes), lwd=2, cex=1.2)
abline(v=eD.out_multi@metadata[["k_means_intercept"]])

hist(log_comb, breaks=200,ylim= c(0,10000), xlab=expression(Log[10]~"(ATAC count)"-"slope * "~Log[10]~"(RNA count)"))
rect(0, 0, eD.out_multi@metadata[["k_means_intercept"]], par("usr")[4], col = "pink", border = NA)
par(new = TRUE)
rect(eD.out_multi@metadata[["k_means_intercept"]], 0 ,12, par("usr")[4], col = "lightblue", border = NA)
par(new = TRUE)
hist(log_comb, breaks=200,ylim= c(0,10000), xlab=expression(Log[10]~"(ATAC count)"-"slope * "~Log[10]~"(RNA count)"))
abline(v=eD.out_multi@metadata[["k_means_intercept"]], col="red")

dev.off()


if (substring(stub, 1,4)=="vale"){
  
  atlas <- read.table( paste0("data/input/", stub, "/FCA_GND1028",  substring(stub, 11,14), "_cell_type.csv"), sep=",", row.names = 1, header=T)
  vale_annotate(atlas, srat_subset, corr_csv, corr_inv_csv)
  
  if (stub=="valentina_8176"   ){

    atlas <- read.table("data/input/valentina_8176/FCA_GND10288176_cell_type_germ.csv", sep=",", row.names = 1, header=T)
    vale_annotate(atlas, srat_subset, corr_csv_fine, corr_inv_csv_fine )
    
  }
  
  
}


srat_subset_res1 <- FindClusters(srat_subset, resolution = 1, random.seed = 0)
pdf(cl_res1)
print(DimPlot(srat_subset_res1, label=T) )
dev.off()


# ATAC clustering
pdf(ATAC_clust)
srat_atac <- srat_subset
DefaultAssay(srat_atac) <- "ATAC"
cl4rna = colnames(srat_subset)[srat_subset$seurat_clusters=="4"]
srat_atac <- RunTFIDF(srat_atac)
srat_atac <- FindTopFeatures(srat_atac, min.cutoff = 'q0')
srat_atac <- RunSVD(srat_atac)
ElbowPlot(srat_atac, n=40)
DepthCor(srat_atac, n=33)
srat_atac <- RunUMAP(object = srat_atac, reduction = 'lsi', dims = 3:30)
srat_atac <- FindNeighbors(object = srat_atac, reduction = 'lsi', dims = 2:30)
srat_atac <- FindClusters(object = srat_atac, verbose = FALSE, algorithm = 3, resolution = 2)
print(DimPlot(object = srat_atac, label = TRUE)+ggtitle("clustering based on atac")) #| DimPlot(object = srat_combined, label = TRUE)+ggtitle("clustering based on rna")  
print(FeaturePlot(srat_atac, features="nCount_RNA")& scale_colour_gradientn(colours = c( "cyan", "deepskyblue", "forestgreen", "darkorange", "darkred"),limits = c(0, 30000)) )
DefaultAssay(srat_atac) <- "SCT"
print(FeaturePlot(srat_atac, order=T, features=c("DCC", "ZP3"), reduction="umap" )& scale_colour_gradientn(colours = c( "cyan", "deepskyblue", "forestgreen", "darkorange", "darkred"))  )
print(FeaturePlot(srat_atac, order=T, features=c("DCC", "ZP3"), reduction="umap" )& scale_colour_gradientn(colours = c( "cyan", "deepskyblue", "forestgreen", "darkorange", "darkred"))  )
print(DimPlot(srat_atac, group.by = "comparison", sizes.highlight=0.1, reduction = "umap", label=T, cols = c("grey", "blue", "red")) )  #+scale_color_manual(labels = c( "eD", "both", "cR",), values = c("blue", "red"))    )


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
write.table(venn_df, venns_per_ATAC_cl, sep = '\t', row.names = T, col.names = T, quote = F)
eD_tpr_only <- unname(unlist(venn_df["% of eD cells",])) - unname(unlist(venn_df["% of common cells",]) )
cR_tpr_only <- unname(unlist(venn_df["% of cR cells",])) - unname(unlist(venn_df["% of common cells",]) )
common <- unname(unlist(venn_df["% of common cells",]))

Values <- matrix(c(cR_tpr_only, common, eD_tpr_only), nrow = 3, ncol = max(as.integer(srat_atac$seurat_clusters)), byrow = TRUE)
barplot(Values, main = "eD vs cR by cluster", names.arg = seq(0, max(as.integer(srat_atac$seurat_clusters))-1) , 
        xlab = "cluster", ylab = "fraction", col = c("grey", "pink", "salmon"))
#legend(25, 0.1, lwd=3, col=c("salmon", "pink", "grey"), lty=c(1,1,1), legend=c("EmptyDrops_multiome ", "common", "cellRanger-arc") )
legend("topleft", bg="white", lwd=3, col=c("salmon", "pink", "grey"), lty=c(1,1,1),
       legend=c("EmptyDrops_multiome ", "common", "cellRanger-arc") )

UMAP_comparison(srat_atac, srat_subset)
UMAP_comparison(srat_subset, srat_atac)


dev.off()
saveRDS(srat_atac, srat_atac_file)


}

