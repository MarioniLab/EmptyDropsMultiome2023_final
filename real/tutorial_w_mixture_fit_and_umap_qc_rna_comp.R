#setwd("./")
library(Seurat)
library(Signac)
library(DropletUtils)
library(ggplot2)
library(ggridges)
library(mixtools)
library(UpSetR)
library(eulerr)
#library(DelayedMatrixStats)
library(EmptyDropsMultiome)
library(dplyr)
#library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
#source("real/vale_annotate.R")
#source("real/clustering_comparison.R")
source("simulations/fcn_for_sim.R")
current_date= paste(unlist(strsplit(as.character(Sys.Date()), "-")), collapse="")
opath <- paste0("data/output/realdata_rna/", current_date)
dir.create(file.path(opath), recursive=TRUE)


ALLFILES <- c(#"pbmc_gran_sorted/raw_feature_bc_matrix.h5",
               "valentina_8176/FCA_GND10288176_raw_feature_bc_matrix.h5",
               "valentina_8177/FCA_GND10288177_raw_feature_bc_matrix.h5"
              # "valentina_8178/FCA_GND10288178_raw_feature_bc_matrix.h5",
              # "valentina_8179/FCA_GND10288179_raw_feature_bc_matrix.h5",
              # "valentina_8180/FCA_GND10288180_raw_feature_bc_matrix.h5"
              
            )

ALL_BARCODES <- c(#"pbmc_gran_sorted/barcodes.tsv.gz",
                  "valentina_8176/barcodes.tsv.gz",
                  "valentina_8177/barcodes.tsv.gz"
                  # "valentina_8178/barcodes.tsv.gz",
                  # "valentina_8179/barcodes.tsv.gz",
                  # "valentina_8180/barcodes.tsv.gz"
)


for (i in seq_along(ALLFILES) ) { 

print("now let's do the next sample")
fname = ALLFILES[i]
sce <- Read10X_h5(file.path("data/input", fname))

stub <- sub("/.*", "", fname, "_qc")
cR_barcodes_file <-  gzfile(file.path("data/input", ALL_BARCODES[i] ))

metadata10x <- read.csv(
  file = paste0("data/input/", stub, "/FCA_GND1028",  substring(stub, 11,14), "_per_barcode_metrics.csv"),
  header = TRUE,
  row.names = 1
)

# create sample specific subdirectory
dir.create(file.path(opath, stub), recursive=TRUE)

# define output file names
eD_multi_tsv <- file.path(opath, stub, paste0(stub, "_eD_multiome.tsv")) 
ffile <- file.path(opath, stub, paste0(stub, "_eD_multiome.pdf"))  
markers_tsv <- file.path(opath, stub, paste0(stub, "_cluster_markers.tsv")) 
eD_metadata <- file.path(opath, stub, paste0("eD_metadata_", stub, ".tsv"))  
venns_per_cl <- file.path(opath, stub, paste0("venns_per_cl_", stub, ".csv"))  
venns_per_cl_pdf <- file.path(opath, stub, paste0("venns_per_cl_", stub, ".pdf"))  
srat_file <- file.path(opath, stub, paste0("srat_", stub, ".rds"))  
print(ffile)
# define count matrices
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

write.table(eD.out_multi,
            paste0(eD_multi_tsv),
            sep = '\t', row.names = T, col.names = T, quote = F)

write.table(eD.out_multi@metadata,
            eD_metadata,
            sep = '\t', row.names = T, col.names = T, quote = F)





cR_cells <- readLines(cR_barcodes_file)

# create srat object
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

DefaultAssay(srat) <- "RNA"

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

union_cells <- union( eD_cells, eD_rna_cells)
eD_minus_rna <- setdiff(eD_rna_cells, eD_cells)
rna_minus_eD <- setdiff( eD_cells, eD_rna_cells )
intersection <- intersect(eD_cells, eD_rna_cells)
srat_subset_rna <- subset(srat, cells=union_cells )
DefaultAssay(srat_subset_rna) <- "RNA"

# add ribo and mito metadata
C <- srat_subset_rna@assays[["RNA"]]
rb.genes <- rownames(C)[grep("^RP[SL]",rownames(C))]
percent.ribo <- colSums(C[rb.genes,])/Matrix::colSums(C)*100
srat_subset_rna <- AddMetaData(srat_subset_rna, percent.ribo, col.name = "percent.ribo")
srat_subset_rna[["percent.mt"]] <- PercentageFeatureSet(srat_subset_rna, pattern = "^MT-")
mito_lim <- median(srat_subset_rna[["percent.mt"]][,1]) + 3* mad(srat_subset_rna[["percent.mt"]][,1])

hist(srat_subset_rna[["TSS.enrichment"]][,1] , breaks=100)

print(VlnPlot(
  object = srat_subset_rna,
  features = c(
    'TSS.enrichment'),
  pt.size = 0.1
)+NoLegend() )


plot.new()
plot.new()
# retained = sum(srat_subset_rna[["percent.ribo"]][,1] < ribo_lim & srat_subset_rna[["percent.mt"]][,1] < mito_lim )
# retained_in_cR = colnames(srat_subset_rna)[srat_subset_rna[["percent.ribo"]][,1] < ribo_lim & srat_subset_rna[["percent.mt"]][,1] < mito_lim & colnames(srat_subset_rna) %in% cR_cells]
# retained_in_eD = colnames(srat_subset_rna)[srat_subset_rna[["percent.ribo"]][,1] < ribo_lim & srat_subset_rna[["percent.mt"]][,1] < mito_lim & colnames(srat_subset) %in% eD_cells]
# rejected = sum(srat_subset_rna[["percent.ribo"]][,1] > ribo_lim | srat_subset_rna[["percent.mt"]][,1] > mito_lim )
# text(x=0.2, y=.1, paste0("retained=", retained))
# text(x=0.2, y=0.2, paste0("rejected=", rejected))
# text(x=0.2, y=0.3, paste0("retained_in_cR=", length(retained_in_cR) ) )
# text(x=0.2, y=0.4, paste0("retained_in_eD=", length(retained_in_eD)) )
# text(x=0.2, y=0.5, paste0("slope=", eD.out_multi@metadata[["k_means_slope"]]))
# text(x=0.2, y=0.6, paste0("intercept=", eD.out_multi@metadata[["k_means_intercept"]] ))

# listInput <- list(retained_eD = retained_in_eD, 
#                   retained_cR = retained_in_cR
# )
# p1 <- upset(fromList(listInput), nsets = 6,, order.by = "freq")
# overlaps <- euler(listInput, shape = "ellipse")
# p2 <- plot(overlaps, 
#            quantities = TRUE,
#            labels = list(font = 4))
# print(p1)
# print(p2)



# subset the srat
srat_subset_rna <- subset(x = srat_subset_rna, subset = percent.mt < mito_lim)
srat_subset_rna <- SCTransform(srat_subset_rna)
srat_subset_rna <- RunPCA(srat_subset_rna, seed.use=42, features = VariableFeatures(object = srat_subset_rna))
print(ElbowPlot(srat_subset_rna, ndims = 50)     )
srat_subset_rna <- FindNeighbors(srat_subset_rna, dims = 1:50)
srat_subset_rna <- FindClusters(srat_subset_rna, resolution = 2, random.seed = 0)
srat_subset_rna <- RunUMAP(srat_subset_rna, dims = 1:50, seed.use = 42)
print(DimPlot(srat_subset_rna, label=T))

srat_subset_rna$comparison <- 1
srat_subset_rna$comparison[ colnames(srat_subset_rna) %in% eD_minus_rna] <- "EmptyDropsMultiome"
srat_subset_rna$comparison[ colnames(srat_subset_rna) %in% intersection] <- "both"
srat_subset_rna$comparison[ colnames(srat_subset_rna) %in% rna_minus_eD] <- "EmptyDrops"

    

retained = length( colnames(srat_subset_rna) )
retained_in_eDmulti = colnames(srat_subset_rna)[ colnames(srat_subset_rna) %in% eD_cells ]
retained_in_eDrna = colnames(srat_subset_rna)[ colnames(srat_subset_rna) %in% eD_rna_cells ]
listInput <- list(eD_after = retained_in_eDmulti, 
                  eDrna_after =retained_in_eDrna
                  )
p1 <- upset(fromList(listInput), nsets = 6,, order.by = "freq")
overlaps <- euler(listInput, shape = "ellipse")
p2 <- plot(overlaps, 
           quantities = TRUE,
           labels = list(font = 4))
print(p1)
print(p2)
    
    
# compute frip_of_excluded & frip_of_cells
FRiP_sub <- srat_subset_rna$atac_peak_region_fragments / srat_subset_rna$atac_fragments * 100
max_frip_of_excluded = max(srat_subset_rna$FRiP[srat_subset_rna$excluded_reason==2]  )
min_frip_of_cells = min(srat_subset_rna$FRiP[srat_subset_rna$is_cell==1]  )

df_FRiP <- data.frame("frip"=srat_subset_rna$FRiP, "cluster"=srat_subset_rna$seurat_clusters)
print(ggplot(df_FRiP, aes(x = frip, y = cluster, height = stat(density))) + 
  geom_density_ridges(stat = "binline", bins = 80, scale = 0.95, draw_baseline = FALSE)+
  theme_ridges(grid = FALSE, center_axis_labels = TRUE)+
  geom_vline(xintercept = max_frip_of_excluded))

  
#saveRDS(srat_subset_rna, paste0(srat_file,"")
saveRDS(srat_subset_rna, srat_file)


# calculate cR/eD venn diagram by cluster
venn_df <- data.frame("0" = c(as.integer(sum(srat_subset_rna$seurat_clusters=="0")),
                              as.integer(sum(srat_subset_rna$seurat_clusters=="0" & colnames(srat_subset_rna) %in% eD_cells)),
                              as.integer(sum(srat_subset_rna$seurat_clusters=="0" & colnames(srat_subset_rna) %in% eD_rna_cells)),
                              sum(srat_subset_rna$seurat_clusters=="0" & colnames(srat_subset_rna) %in% eD_rna_cells & colnames(srat_subset_rna) %in% eD_cells),
                              sum(srat_subset_rna$seurat_clusters=="0" & colnames(srat_subset_rna) %in% eD_cells)/sum(srat_subset_rna$seurat_clusters=="0"),
                              sum(srat_subset_rna$seurat_clusters=="0" & colnames(srat_subset_rna) %in% eD_rna_cells)/sum(srat_subset_rna$seurat_clusters=="0"),
                              sum(srat_subset_rna$seurat_clusters=="0" & colnames(srat_subset_rna) %in% eD_rna_cells & colnames(srat_subset_rna) %in% eD_cells)/sum(srat_subset_rna$seurat_clusters=="0")
                                  ))
for ( cl in c(1:(max(as.integer(srat_subset_rna$seurat_clusters))-1))   ){
  char_cl = as.character(cl)
  temp_df <- data.frame("new" = c(sum(srat_subset_rna$seurat_clusters==char_cl),
                                           sum(srat_subset_rna$seurat_clusters==char_cl & colnames(srat_subset_rna) %in% eD_cells),
                                           sum(srat_subset_rna$seurat_clusters==char_cl & colnames(srat_subset_rna) %in% eD_rna_cells),
                                           sum(srat_subset_rna$seurat_clusters==char_cl & colnames(srat_subset_rna) %in% eD_rna_cells & colnames(srat_subset_rna) %in% eD_cells),
                                           sum(srat_subset_rna$seurat_clusters==char_cl & colnames(srat_subset_rna) %in% eD_cells)/sum(srat_subset_rna$seurat_clusters==char_cl),
                                           sum(srat_subset_rna$seurat_clusters==char_cl & colnames(srat_subset_rna) %in% eD_rna_cells)/sum(srat_subset_rna$seurat_clusters==char_cl),
                                           sum(srat_subset_rna$seurat_clusters==char_cl & colnames(srat_subset_rna) %in% eD_rna_cells & colnames(srat_subset_rna) %in% eD_cells)/sum(srat_subset_rna$seurat_clusters==char_cl)
  ))
  venn_df <- cbind(venn_df, temp_df)
}
colnames(venn_df) <- as.character(c(0:(max(as.integer(srat_subset_rna$seurat_clusters))-1) ))
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
Values <- matrix(c(cR_tpr_only, common, eD_tpr_only), nrow = 3, ncol = max(as.integer(srat_subset_rna$seurat_clusters)), byrow = TRUE)
barplot(Values, main = "eD vs cR by cluster", names.arg = seq(0, max(as.integer(srat_subset_rna$seurat_clusters))-1) , 
        xlab = "cluster", ylab = "fraction", col = c("grey", "pink", "salmon"))
#legend(25, 0.1, lwd=3, col=c("salmon", "pink", "grey"), lty=c(1,1,1), legend=c("EmptyDrops_multiome ", "common", "cellRanger-arc") )
legend("topleft", bg="white", lwd=3, col=c("salmon", "pink", "grey"), lty=c(1,1,1),
       legend=c("EmptyDrops_multiome ", "common", "cellRanger-arc") ) 
    

dev.off()





}

