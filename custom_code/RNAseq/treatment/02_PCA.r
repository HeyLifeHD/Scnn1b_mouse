#libraries
library(DESeq2)
library(ChIPseeker)
require(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
library(dplyr)
library(pheatmap)
library(ggpubr)
library(limma)

#folder
base.dir <- "/home/heyj/c010-datasets/Internal/COPD/NaturComm/RNAseq/treatment"
base_results.dir <- file.path(base.dir, "results")
results.dir<- file.path(base_results.dir , "tables")
PreDE.dir <- file.path(base_results.dir,"PreDE")
PostDE.dir <- file.path(base_results.dir,"PostDE")

#Data Read-in
#DESeq2 Analysis
dds <- readRDS(file = file.path(results.dir,"dds.rds"))
anno <- colData(dds)
#Extracting transformed values 
rld_RNA <- rlog(dds, blind=FALSE)
saveRDS(rld_RNA,file = file.path(results.dir,"rld_all.rds"))
rlog_counts <- assay(rld)
#correct for batch effects
BR<-removeBatchEffect(rlog_counts, batch=colData(rld)$Date)
rld_b <- SummarizedExperiment(assays =BR , colData=colData(rld))
rld_b <- DESeqTransform(rld_b)
saveRDS(rld_b,file =file.path(results.dir, "rld_b_date.rds"))

#pca analysis
pca <- prcomp(t(assay(rld_b)),
                 center = TRUE,
                 scale. = F) 

summary(pca)
x <- pca$x
pheno <- colData(rld_b_Date_RNA)
x<- as.data.frame(cbind(x , pheno))

pdf(file.path("c010-datasets/Internal/COPD/LPS_Integration/PCA", "PCA_RNA_rld_b_date.pdf"), height=3.5, width=5)
ggscatter(x, x="PC1", y="PC2",
          color = "Genotype", shape = "Treatment", size=3,
          ellipse = F , mean.point = FALSE, palette= c("#EFC000FF","#868686FF"),
          star.plot = F, xlab=(paste0("PC1: ", round(summary(pca)$importance[2,1]*100,2), "% variance")), ylab=(paste0("PC2: ", round(summary(pca)$importance[2,2]*100,2), "% variance")))+theme(legend.position="right")
dev.off()
