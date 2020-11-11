
#Directories
base.dir <- "/home/heyj/c010-datasets/Internal/COPD/NaturComm/ATAC/treatment"
analysis.dir <- file.path(base.dir, "analysis")
data.dir <- "icgc/dkfzlsdf/analysis/C010/betaENaC/ATAC/LPS_experiment/processing/"

#load data
sample_anno<- readRDS(file.path(analysis.dir,"data", "sample_anno.rds"))
counts_raw<- readRDS(file.path(analysis.dir, "allPeaks_raw.rds"))

#with interactive term
dds <- DESeqDataSetFromMatrix(countData = mcols(counts_raw), 
                              colData = sample_anno, 
                              design = ~  Replicate + Factor, rowRanges=counts_raw)
#for Genotype comparison
dds <- estimateSizeFactors(dds)
dim(dds)
#Running the differential expression 
dds <- DESeq(dds)
saveRDS(dds,file = file.path(analysis.dir,"dds_group.rds"))

#extract transformed values
rld_ATAC <- rlog(dds, blind=FALSE)
saveRDS(rld_ATAC,file = file.path(analysis.dir,"rld.rds"))
rld_ATAC<- readRDS(file.path(analysis.dir,"rld.rds"))

#pca analysis
rld_ATAC_df <- as.data.frame(assay(rld_ATAC))
pca <- prcomp(t(rld_ATAC_df),
                 center = TRUE,
                 scale. = F) 

summary(pca)
x <- pca$x
pheno <- colData(rld_ATAC)
colnames(pheno)<- c("SampleID" ,  "Genotype"   ,  "Treatment"  ,"Replicate"  ,"Date"     , "factor"   ,  "Peaks"     , "bamReads"  , "bams"   ,    "peaks"     ,"Factor"  ,   "sizeFactor")
x<- as.data.frame(cbind(x , pheno))

pdf(file.path(base.dir, "PCA_ATAC.pdf"), height=3.5, width=5)
ggscatter(x, x="PC1", y="PC2",
          color = "Genotype", shape = "Treatment", size=3,
          ellipse = F , mean.point = FALSE, palette= c("#EFC000FF","#868686FF"),
          star.plot = F, xlab=(paste0("PC1: ", round(summary(pca)$importance[2,1]*100,2), "% variance")), ylab=(paste0("PC2: ", round(summary(pca)$importance[2,2]*100,2), "% variance")))+theme(legend.position="right")
dev.off()
