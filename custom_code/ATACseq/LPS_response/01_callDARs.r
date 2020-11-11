
#Directories
#NGS-machine
base.dir <- "/home/heyj/c010-datasets/Internal/COPD/NaturComm/ATAC/treatment"
data.dir <- "/icgc/dkfzlsdf/analysis/C010/betaENaC/ATAC/BL_experiment/processing/"
analysis.dir <- file.path(base.dir, "analysis")
dir.create(analysis.dir, recursive=TRUE)



#libraries
require(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
library(ChIPpeakAnno)
library(reshape2)
library(ggpubr)
library(DESeq2)
#load data
sample_anno<- readRDS(file.path(analysis.dir,"data", "sample_anno.rds"))
counts_raw<- readRDS(file.path(analysis.dir, "allPeaks_raw.rds"))
colnames(sample_anno)[2]<-"Genotype"

#create dds
model.matrix(~ Date + Genotype + Treatment + Genotype:Treatment , sample_anno)
dds <- DESeqDataSetFromMatrix(countData = mcols(counts_raw), rowRanges=counts_raw,
                              colData = sample_anno, 
                              design = ~ Date + Genotype + Treatment + Genotype:Treatment )
#Estimating size factors
dds <- estimateSizeFactors(dds)
#Running the differential expression 
dds <- DESeq(dds)
saveRDS(dds,file = file.path(analysis.dir, "dds.rds"))


#Set specifications
alpha <- 0.05 #set FDR cutoff
lfc <- 0 #set logfold2 cutoff

#run the differential accessibility analysis 
results_mediumtg_vs_LPStg_inter<- results(dds, contrast=c("Treatment","medium", "LPS"), alpha = alpha, lfcThreshold = lfc, format="GRanges") #extract results of wanted comparison
results_mediumwt_vs_LPSwt_inter<- results(dds,list(c("Treatment_medium_vs_LPS", "Genotypewt.Treatmentmedium")) , alpha = alpha, lfcThreshold = lfc, format="GRanges") #extract results of wanted comparison
results_Interaction <-results(dds, name="Genotypewt.Treatmentmedium", alpha = alpha, lfcThreshold = lfc, format="GRanges")
results_Genotype <-results(dds, name="Genotype_wt_vs_tg", alpha = alpha, lfcThreshold = lfc, format="GRanges")
results_treatment <- results(dds, name="Treatment_medium_vs_LPS", alpha = alpha, lfcThreshold = lfc, format="GRanges")
#summmaries
summary(results_mediumtg_vs_LPStg_inter) 
summary(results_mediumwt_vs_LPSwt_inter)
summary(results_Interaction) 
summary(results_Genotype) 
summary(results_treatment)

#create list
results <- list(results_mediumtg_vs_LPStg_inter, results_mediumwt_vs_LPSwt_inter, results_Interaction,results_Genotype, results_treatment )
names(results)<- c("tg_medium_vs_tg_LPS", "wt_medium_vs_wt_LPS", "Interaction","Genotype", "Treatment")


#annotate list
DE_list <- list(NULL)
for (i in names(results)){
    #extract DARs
    temp <- results[[i]]
    if(is.null(temp)) {}
    if(length(temp)==1) {
    temp_anno<- annotatePeak(peak= temp, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Mm.eg.db")
    #make df
    temp_anno_df <- as.data.frame(temp_anno)
    temp_anno_gr <- makeGRangesFromDataFrame(temp_anno_df, keep.extra.columns=TRUE)
    #make a list
    DE_list[[i]]<- temp_anno_gr
    }
    if(length(temp)>1){
        #annotate DE regions
    temp_anno<- annotatePeak(peak= temp, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Mm.eg.db")
    #make df
    temp_anno_df <- as.data.frame(temp_anno)
    temp_anno_gr <- makeGRangesFromDataFrame(temp_anno_df, keep.extra.columns=TRUE)
    #make a list
    DE_list[[i]]<- temp_anno_gr
    }
}

DE_list<- lapply(DE_list, function(x){
    x <- x[order(x$padj, decreasing=FALSE),]
    x
})
saveRDS(DE_list, file.path(analysis.dir, "de_list_complete.rds"))