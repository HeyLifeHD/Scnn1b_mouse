#libraries
library(DESeq2)
library(ggpubr)
library(rafalib)
library(pheatmap)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(biomaRt)
library(knitr)
library(clusterProfiler)
library(Glimma)
library(limma)
library(edgeR)

#directories
base.dir <- "/home/heyj/c010-datasets/Internal/COPD/NaturComm/RNAseq/treatment"
data.dir <- file.path("/home/heyj/c010-datasets/Internal/COPD/RNASeq/190612_3rep_MedLPS", "data")
base_results.dir <- file.path(base.dir, "results")
results.dir<- file.path(base_results.dir , "tables")
PreDE.dir <- file.path(base_results.dir,"PreDE")
PostDE.dir <- file.path(base_results.dir,"PostDE")

#Read in Data
dds <- readRDS(file = file.path(results.dir,"dds.rds"))
anno <- colData(dds)

#Take a look at design and annotation of samples
design(dds)

#extract with results
#Set specifications
alpha <- 0.1 #set FDR cutoff
lfc <- 0 #set logfold2 cutoff

#Running the differential expression 
#set up all possible contrasts
contrasts <- as.data.frame(combn(as.character(unique(colData(dds)$group)), 2))
results <- list()
for (i in 1:length(contrasts)) {
  results[[i]]<- results(dds, contrast=c("group",as.character(contrasts[1,i]),  as.character(contrasts[2,i])), alpha = alpha, lfcThreshold = lfc) #extract results of wanted comparison
  print(paste0( as.character(contrasts[1,i]), "_vs_",  as.character(contrasts[2,i])))
  print(summary(results[[i]]))
}

nam<-vector()
for (i in 1:length(contrasts)) {
nam[i]<- paste0( as.character(contrasts[1,i]), "_vs_",  as.character(contrasts[2,i]))
}
names(results)<- nam

#annotate samples
#Make a dataframe out of it
DEG_results_list <- lapply(results, function(x){
  x<- as.data.frame(x)
  x
})
#create ensembl annotation
DEG_results_list <- lapply(DEG_results_list, function(x){
  x$ensembl <- sapply(strsplit(rownames(x) ,".", fixed=TRUE),`[`, 1)
  x
})
#add symbol annotation
DEG_results_list <- lapply(DEG_results_list, function(x){
  x$symbol<- mapIds(org.Mm.eg.db, keys=x$ensembl, keytype ="ENSEMBL", column = "SYMBOL", multiVals = "first" )
  x
})
#add entrezgene annotation
DEG_results_list <- lapply(DEG_results_list, function(x){
  x$entrezgene<- mapIds(org.Mm.eg.db, keys=x$ensembl, keytype ="ENSEMBL", column = "ENTREZID", multiVals = "first" )
  x
})
#add gene name annotation
DEG_results_list <- lapply(DEG_results_list, function(x){
  x$genename<- mapIds(org.Mm.eg.db, keys=x$ensembl, keytype ="ENSEMBL", column = "GENENAME", multiVals = "first" )
  x
})
#order by padjusted value
DEG_results_list <- lapply(DEG_results_list, function(x){
  x <- x[order(x$padj),]
  x
})
names(DEG_results_list)<-names(results)

#save lists
saveRDS(DEG_results_list, file.path(PostDE.dir, "DEG_results_group_list.rds"))
