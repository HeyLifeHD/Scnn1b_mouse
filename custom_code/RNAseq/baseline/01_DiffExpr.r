
#libraries
library(DESeq2)
library(biomaRt)
library(ggplot2)
library(limma)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(IHW)
library(biomaRt)

#folders
##output
base.dir <- "/home/heyj/c010-datasets/Internal/COPD/NaturComm/"
base_results.dir <- file.path(base.dir, "03_results")
results.dir<- file.path(base_results.dir , "tables")
PreDE.dir <- file.path(base_results.dir,"PreDE")
PostDE.dir <- file.path(base_results.dir,"PostDE")
##input
data.dir <- file.path("/Users/c010/Documents/UNI/Phd/COPD Project/RNA_Seq/180503_AM/", "01_data")

#load data
CT <- readRDS(file.path(data.dir, "sample_annotation_new.rds"))
GC <-readRDS(file.path(data.dir, "gene_count_new.rds"))
dds <- DESeqDataSetFromMatrix(countData = GC, 
                              colData = CT, design = ~ batch+ Genotype)

#Estimating size factors
dds <- estimateSizeFactors(dds)

#Filter genes which are only expressed in 1 sample
idx <- rowSums( counts(dds, normalized=TRUE) >= 1 ) >= 3
dds <- dds[idx,]
dim(dds)

#Estimate dd
#for Genotype comparison
dds <- estimateDispersions(dds)

pdf(file.path(PreDE.dir, "Dispersion(>1in2).pdf"))
plotDispEsts(dds, main="Dispersion plot")
dev.off()

#Running the differential expression 
#for Genotype comparison
dds <- DESeq(dds,minReplicatesForReplace=3)
dds<-replaceOutliers(dds,minReplicates=3)
saveRDS(dds,file = file.path(results.dir,"dds.rds"))

#diff expression analysis
result <- results(dds, contrast = c("Genotype",  "tg",  "wt"), alpha=0.1, lfcThreshold = 0)
summary(result)
saveRDS(result,file = file.path(results.dir,"180404_result.rds"))
resultdf<- as.data.frame(result)
resultdf$DE <- NA
resultdf$DE <- ifelse(abs(resultdf$log2FoldChange)>.5 &resultdf$padj  <0.1, "Sig.DE","notDE" )
table(resultdf$DE)
#shrink log fold changes
resultshrink <- lfcShrink(dds,contrast = c("Genotype",  "tg",  "wt"), res=result )
summary(resultshrink)
saveRDS(resultshrink,file = file.path(results.dir,"180404_result_shrink.rds"))
