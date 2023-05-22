setwd("/mnt/beegfs6/home3/ahringer/em613/analysis/multiomics/emptryDrops_multiome2023")
library(eDv2)
library(DropletUtils)
library(Matrix)
library(rstudioapi)
library(Seurat)
library(mixtools)
library(tidyverse)
source("simulations/fcn_for_sim_more_empties.R")
library("optparse")



opath <- "data/output/sim_more_empties/results-sim"
ppath <- "data/output/sim_more_empties/pics-sim"

dir.create(opath,recursive=TRUE)
dir.create(ppath,recursive=TRUE)

set.seed(73953024)
#set.seed(10)


ALLFILES <- c("pbmc_gran_sorted/raw_feature_bc_matrix.h5",
              #"gastr_d4_A/raw_feature_bc_matrix.h5",
              #"gastr_d4_B/gastr_d4_multiome_B_L001_raw_feature_bc_matrix.h5",
              "valentina_FCA_GND10288180/raw_feature_bc_matrix" #,
              #"valentina_FCA_GND10288177/FCA_GND10288177_raw_feature_bc_matrix.h5"
)

# choose specific file
option_list = list(
  make_option(c("-n", "--number"), type="character", default=NULL, 
              help="which sample", metavar="character")
); 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
# define inputs/outputs and useful functions
i <- opt$number 
ALLFILES <- c( ALLFILES[i] )


for (fname in ALLFILES) {
  #fname =  ALLFILES[i]
  #sce <- Read10X_h5(file.path("../../..", "data/input", fname))
  #sce <- Read10X_h5(fname)
  stub <- sub("/.*", "", fname)
  if (stub=="valentina_FCA_GND10288180"){
    sce <- Read10X(file.path("data/input", fname))
  }else{
    sce <- Read10X_h5(file.path("data/input", fname))
  }
  
  ofile <- file.path(opath, paste0(stub, "eD-multiome-res.tsv"))
  ffile <- file.path(opath, paste0(stub, "-fdr.tsv"))
  unlinker <- TRUE 
  
  raw.mat <- sce[["Gene Expression"]]
  raw.mat_atac <- sce[["Peaks"]]
  
  for (g1 in c(1000, 5000)){ 
    for (g2 in c(1000, 5000)) {
      first <- TRUE 
      print(g1)
      print(g2)
      
      for (it in 1:1) { 
        print(it)
        #out <- SIMFUN_multiome(raw.mat, raw.mat_atac, group1=g1, group2=g2, down.rate=0.02)
        # for gastr and pbmc down.rate=0.02 (for atac; and twice this for rna)
        #out <- SIMFUN_multiome_v2(raw.mat, raw.mat_atac, group1=g1, group2=g2, down.rate=0.02)
        if (stub=="valentina_FCA_GND10288180"){
          print("correct")
          #out <- SIMFUN_multiome_v2(raw.mat, raw.mat_atac, group1=g1, group2=g2, reorder.rate=0.001, down.rate.rna=0.2, down.rate.atac=0.2)  #0.1
          #out <- SIMFUN_multiome_v3(raw.mat, raw.mat_atac, group1=g1, group2=g2, reorder.rate=0.001, down.rate.rna=0.2, down.rate.atac=0.2, lower_atac=60, lower_rna=110)  
          #out <- SIMFUN_multiome_no_scrambling(raw.mat, raw.mat_atac, group1=g1, group2=g2, reorder.rate=0.1, down.rate.rna=0.7, down.rate.atac=0.7)  
          out <- SIMFUN_multiome_no_scrambling_vale_more_empties(raw.mat, raw.mat_atac, group1=g1, group2=g2, reorder.rate=0.1, down.rate.rna=0.05, down.rate.atac=0.05)  
        }else if(stub=="valentina_FCA_GND10288177"){
          out <- SIMFUN_multiome_no_scrambling(raw.mat, raw.mat_atac, group1=g1, group2=g2, reorder.rate=0.1, down.rate.rna=0.1, down.rate.atac=0.1)  
        }else{
          #out <- SIMFUN_multiome_v2(raw.mat, raw.mat_atac, group1=g1, group2=g2, down.rate.rna=0.04, down.rate.atac=0.02)
          
          #out <- SIMFUN_multiome_no_scrambling(raw.mat, raw.mat_atac, group1=g1, group2=g2, reorder.rate=0.1, down.rate.rna=0.04, down.rate.atac=0.02)  
          out <- SIMFUN_multiome_no_scrambling_more_empties(raw.mat, raw.mat_atac, group1=g1, group2=g2, reorder.rate=0.1, down.rate.rna=0.04, down.rate.atac=0.02)  
        }
        dev.off()
        
        final_rna <- out$counts_rna
        final_atac <- out$counts_atac
        totals_rna <- colSums(final_rna)
        totals_atac <- colSums(final_atac)
        
        plot(log10( colSums(final_atac) + 0.1), log10( colSums(final_rna)  + 0.1),
             #     xlim=c(0,1510) , ylim=c(0,8),
             pch=".",
             # cex=2,
             col=factor(out$identity ),
             xlab="log10(atac_count)", ylab="log10(rna_count)",
             main=paste0("ground truth")
        )+theme(axis.title.x=element_text(size=14,face="bold"),
                axis.title.y=element_text(size=14,face="bold"))
        
        if (first) {
          o <- order(totals_rna, decreasing=TRUE)
          r <- seq_along(o)
          ot <- totals_rna[o]
          oi <- out$identity[o]
          
          #pdf(file.path("../../..", ppath, paste0(stub, "_", g1, "_", g2, "-sim.pdf")), width=10, height=8)
          pdf(file.path( ppath, paste0(stub, "_multiomics_", g1, "_", g2, "-sim.pdf")), width=10, height=8)
          par(mfrow=c(2,2), cex.axis=1.2, cex.lab=1.4, cex.main=1.4, mar=c(4.1, 4.1, 2.1, 1.1))
          plotBarcodes(r, ot, pch=16, main="All barcodes")
          plotBarcodes(r, ot, pch=16, subset=(oi==0), col="grey80", main="Empty droplets")
          plotBarcodes(r, ot, pch=16, subset=(oi==1), col="blue", main="Large cells")
          plotBarcodes(r, ot, pch=16, subset=(oi==2), col="red", main="Small cells")
          dev.off()
          
          
          
        }
        
        # Testing emptyDrops.
        if (stub=="valentina_FCA_GND10288180"){
          lower_rna = NULL
          barhop_rna = NULL
          lower_atac = 60
          barhop_atac = 3
        }else if(stub=="valentina_FCA_GND10288177"){
          lower_rna = 100
          barhop_rna = 10
          lower_atac = 90
          barhop_atac = 20          
        }else if(stub=="pbmc_gran_sorted" ){
          lower_rna = 150
          barhop_rna = 10
          lower_atac = 110
          barhop_atac = 20          
        }else{
          lower_rna = NULL
          barhop_rna = NULL
          lower_atac = NULL
          barhop_atac = NULL        
        }
        
        # # give index names to missing colnames
        start_time <- Sys.time()
        
        pdf(file.path( ppath, paste0(stub, "_eD_multiome_", g1, "_", g2, "-sim.pdf")), width=10, height=8)
        eD.out_multi <- emptydrops_multiome(final_rna, lower_rna, barhop_rna, final_atac, lower_atac, barhop_atac, niter_rna = 20000, niter_atac = 20000 )  #niter_rna = 10000, niter_atac = 20000
        
        par(mfrow=c(1,1))
        
        end_time <- Sys.time()
        print(end_time - start_time)
        
        
        is.sig <- eD.out_multi$FDR_multi <= 0.001 & !is.na(eD.out_multi$FDR_multi)
        
        # rownames of emptydrops output are shuffled compared to colnames of the simulation. Shuffle them back.
        identities = out$identity[rownames(eD.out_multi)]
        emp.res <- assessMethod(is.sig, identities)
        
        plot(log10(eD.out_multi$Total_chromatin + 0.1), log10(eD.out_multi$Total_RNA + 0.1),
             #     xlim=c(0,1510) , ylim=c(0,8),
             pch=".",
             # cex=2,
             col=factor(identities ),
             xlab="log10(atac_count)", ylab="log10(rna_count)",
             main=paste0("ground truth")
        )+theme(axis.title.x=element_text(size=14,face="bold"),
                axis.title.y=element_text(size=14,face="bold"))
        
        
        plot(log10(eD.out_multi$Total_chromatin + 0.1), log10(eD.out_multi$Total_RNA + 0.1),
             #     xlim=c(0,1510) , ylim=c(0,8),
             pch=".",
             # cex=2,
             col=factor(eD.out_multi$FDR_RNA <= 0.001 & !is.na(eD.out_multi$FDR_RNA) ),
             xlab="log10(atac_count)", ylab="log10(rna_count)",
             main=paste0("eDrna")
        )+theme(axis.title.x=element_text(size=14,face="bold"),
                axis.title.y=element_text(size=14,face="bold"))
        abline(a=eD.out_multi@metadata["k_means_intercept"][[1]],b=eD.out_multi@metadata["k_means_slope"][[1]], col="blue")
        # abline(a=equation_parallel[1],b=equation_parallel[2], col="blue")
        
        dev.off()
        
        write.table(eD.out_multi,  file.path( ppath, paste0(stub, "_eD_multiome_", g1, "_", g2, "-output.csv")), 
                    sep=",", row.names = T, col.names = T, quote = F  )
        write.table(eD.out_multi@metadata,  file.path( ppath, paste0(stub, "_eD_multiome_", g1, "_", g2, "-metadata.csv")), 
                    sep=",", row.names = T, col.names = T, quote = F  )
        
        
        # Using the CellRanger-arc approach.
        cRarc.res <- assessMethod(eD.out_multi$k_means, identities)
        
        # cell calling based on eD-rna.
        eDrna.res <- assessMethod(eD.out_multi$FDR_RNA <= 0.001 & !is.na(eD.out_multi$FDR_RNA), identities)
        
        # cell calling based on eD-atac.
        eDatac.res <- assessMethod(eD.out_multi$FDR_chromatin <= 0.001 & !is.na(eD.out_multi$FDR_chromatin), identities)
        
        
        # Saving the default assessment results to file.
        write.table(cbind(G1Size=g1, G2Size=g2,
                          Method=c("eD-multiome", "eD-rna", "eD-atac", "CellRanger-arc"),
                          rbind(emp.res, eDrna.res, eDatac.res, cRarc.res)),
                    #file=file.path("../../..",ofile), append=!unlinker, col.names=unlinker,
                    file=file.path(ofile), append=!unlinker, col.names=unlinker,
                    row.names=FALSE, quote=FALSE, sep="\t")
        
        if (first) {
          # Creating a ROC plot.
          #hard_cut <- log(eD.out_multi$Total_RNA + 1 + eD.out_multi$Total_chromatin)
          #hard_cut <- log(eD.out_multi$Total_RNA + 1 - eD.out_multi@metadata["k_means_slope"][[1]] * eD.out_multi$Total_chromatin)
          hard_cut <- log10(eD.out_multi$Total_RNA + 1) - eD.out_multi@metadata["k_means_slope"][[1]] * log10(eD.out_multi$Total_chromatin+1)
          
          emp.roc <- createRocPts(log(eD.out_multi$FDR_multi+1e-3), identities)
          lib.roc <- createRocPts(-hard_cut, identities)
          # eDrna.roc <- createRocPts(log(eD.out_multi$FDR_RNA+1e-3), identities)
          emp.fdr <- emp.roc[,"0"]/rowSums(emp.roc)
          emp.tp1 <- emp.roc[,"1"]/g1
          emp.tp2 <- emp.roc[,"2"]/g2
          lib.fdr <- lib.roc[,"0"]/rowSums(lib.roc)
          lib.tp1 <- lib.roc[,"1"]/g1
          lib.tp2 <- lib.roc[,"2"]/g2
          
          emp.col <- colors["eD-multiome"]
          lib.col <- "grey50"
          
          pdf(file.path(ppath, paste0(stub, "_eD_multiome_", "_", g1, "_", g2, "-roc.pdf")), width=8, height=6)
          par(mar=c(5.1, 4.1, 4.1, 12.1))
          plot(c(0, emp.fdr), c(0, emp.tp1), xlim=c(0, 0.07), ylim=c(0, 1), type="l", col=emp.col, lwd=3, main=stub,
               xlab="Observed FDR", ylab="Recall", cex.axis=1.2, cex.lab=1.4, cex.main=1.4)
          lines(c(0, lib.fdr), c(0, lib.tp1), col=lib.col,lwd=3)
          lines(c(0, emp.fdr), c(0, emp.tp2), col=emp.col, lwd=3, lty=2)
          lines(c(0, lib.fdr), c(0, lib.tp2), col=lib.col, lwd=3, lty=2)
          
          par(xpd=TRUE)
          legend(0.035, 1, lwd=3, col=rep(c(emp.col, lib.col), each=2), lty=rep(1:2, 2),
                 legend=c("EmptyDrops (large)", "EmptyDrops (small)", "Total count (large)", "Total count (small)"))
          dev.off()
        }
        
        # Reporting on FDR control.
        all.thresholds <- c(0.001, 0.002, 0.005, 0.01) 
        collected <- numeric(length(all.thresholds))
        for (i in seq_along(all.thresholds)) {
          collected[[i]] <- assessMethod(  eD.out_multi$FDR_multi<= all.thresholds[i], identities )["FDR"]
        }
        
        collected <- rbind(collected)
        colnames(collected) <- all.thresholds
        write.table(cbind(G1Size=g1, G2Size=g2, collected), file=ffile, append=!unlinker,
                    col.names=unlinker, row.names=FALSE, quote=FALSE, sep="\t")
        
        unlinker <- FALSE
        first <- FALSE
      }
    }
  }
}

