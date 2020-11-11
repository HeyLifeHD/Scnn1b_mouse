
#Directories
#PC-09
base.dir <- "/home/heyj/c010-datasets/Internal/COPD/NaturComm/ATAC/"
data.dir <- "/icgc/dkfzlsdf/analysis/C010/betaENaC/ATAC/BL_experiment/processing/"
analysis.dir <- file.path(base.dir, "analysis")
dir.create(analysis.dir, recursive=TRUE)

#Libraries
library(dplyr)
library(DiffBind)
library(ChIPseeker)
library(rtracklayer)
library(pheatmap)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

#load in sample anno
sample_anno <- data.frame(SampleID=c("AM_wt_rep1","AM_wt_rep2", "AM_wt_rep3", 
        "AM_tg_rep1", "AM_tg_rep2", "AM_tg_rep3", "AM_tg_rep4"), 
        Treatment=c(rep("none",3),rep("none",4)),
        Tissue=c( rep("wt",3),rep("tg",4)),
        Factor=c( rep("wt_baseline",3),rep("tg_baseline",4)),
        Replicate=c(1,2,4,1,2,3,4), lane=c(1,1,3,1,1,2,3),
        bamReads=NA,
        Peaks=NA
    )
rownames(sample_anno)<-sample_anno$SampleID

#find bam and bed files
#bams
bams_BL <- dir(path = file.path(data.dir), full.names = TRUE, recursive=TRUE,
             pattern = ".merged.nodup.bam")
bams_BL <- bams_BL[grep("bai", bams_BL, invert=TRUE)]
bams_BL <- bams_BL[grep("glob", bams_BL, invert=TRUE)]
bams_BL <- bams_BL[grep("input", bams_BL, invert=TRUE)]
bams_BL <- c(bams_BL[5:7], bams_BL[1:4])
bams <- bams_BL
#beds
beds_BL <- dir(path = file.path(data.dir), full.names = TRUE, recursive=TRUE,
             pattern = ".merged.nodup.tn5.pval0.01.300K.bfilt.narrowPeak.gz")
beds_BL <- beds_BL[grep("glob", beds_BL, invert=TRUE)]
beds_BL <- c(beds_BL[5:7], beds_BL[1:4])
beds <- c(beds_BL)
#join everything in sample_anno
sample_anno$bamReads <- bams
sample_anno$Peaks <- beds

#saving
dir.create(file.path(analysis.dir,"data"),recursive=TRUE)
saveRDS(sample_anno, file.path(analysis.dir,"data", "sample_anno.rds"))
write.table(sample_anno, file.path(analysis.dir,"data", "sample_anno.tsv"), quote=FALSE, sep="\t", row.names=FALSE)

#create dataset
dataset <-dba(sampleSheet=sample_anno, peakCaller="macs")

#count reads for each peak
dataset<- dba.count(dataset,  bParallel=T)
saveRDS(dataset, file.path(analysis.dir,"Dataset_count.rds"))

#set up the contrast DAR analysis
dataset<- dba.contrast(dataset, categories=DBA_TISSUE, block=DBA_REPLICATE, minMembers=2)

#run DAR analysis
dataset <- dba.analyze(dataset,method=DBA_ALL_METHODS, bParallel=F)
saveRDS(dataset, file.path(analysis.dir, "diffbind.rds"))

#get complete peak set TPM
counts <- dba.peakset(dataset, bRetrieve=TRUE)
saveRDS(counts, file.path(analysis.dir, "countsperpeak.rds"))
write.table(as.data.frame(counts), file.path(analysis.dir, "allPeaks.txt"),quote=F)
export.bed(counts, file.path(analysis.dir, "allPeaks.bed"))

#get complete peak set raw counts 
dataset_rawCounts <-  dba.count(dataset,peaks=NULL, score=DBA_SCORE_READS)
counts_raw <- dba.peakset(dataset_rawCounts, bRetrieve=TRUE,score=DBA_SCORE_READS)
write.table(as.data.frame(counts_raw), file.path(analysis.dir, "allPeaks_raw.txt"),quote=F)
export.bed(counts_raw, file.path(analysis.dir, "allPeaks_raw.bed"))
saveRDS(counts_raw, file.path(analysis.dir, "allPeaks_raw.rds"))

#extact DAR regions for each of the comparisons
contrasts <- data.frame(contrasts=1:length(unlist(lapply(dataset$contrasts, function(x)x$name1))), group1=c(unlist(lapply(dataset$contrasts, function(x)x$name1))),
group2= c(unlist(lapply(dataset$contrasts, function(x)x$name2))))
contrasts$comparison <- paste0(contrasts$group1, "_vs_", contrasts$group2)
#create folders
lapply(c("Anno", "MAs", "Volcanos", "Heatmaps", "Boxplots"), function(x)dir.create(file.path(analysis.dir,x)))

#write for loop to extract all possible comparisons in different plots
DE_list <- list(NULL)
for (i in 1:nrow(contrasts)){
    #extract DARs
    #specific for this mixed analysis
    if(i <=2){
        temp <- dba.report(dataset, contrast=i ,method =DBA_EDGER_BLOCK, th=0.05)
    } else{
        temp <- dba.report(dataset, contrast=i ,method =DBA_EDGER_BLOCK, th=0.05)
    }
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
names(DE_list) <- contrasts$comparison
saveRDS(DE_list, file.path(analysis.dir, "de_list.rds"))
DE_list <- readRDS(file.path(analysis.dir, "de_list.rds"))

#retrieve complete list
DE_list_complete <- list(NULL)
for (i in 1:nrow(contrasts)){
    #specific for this mixed analysis
    if(i <=2){
        temp <- dba.report(dataset, contrast=i ,method =DBA_EDGER_BLOCK, th=1)
    } else{
        temp <- dba.report(dataset, contrast=i ,method =DBA_EDGER_BLOCK, th=1)
    }    
    temp_anno<- annotatePeak(peak= temp, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Mm.eg.db")
    #make df
    temp_anno_df <- as.data.frame(temp_anno)
    temp_anno_gr <- makeGRangesFromDataFrame(temp_anno_df, keep.extra.columns=TRUE)
    #make a list
    DE_list_complete[[i]]<- temp_anno_gr
    }
#order and join with counts
counts <- sortSeqlevels(counts)
counts <- sort(counts)
DE_list_complete <- lapply(DE_list_complete, function(x){
    x <- sortSeqlevels(x)
    x <- sort(x)
    mcols(x)<- cbind(mcols(x), mcols(counts))
    x
    x <- x[order(x$FDR, decreasing=F),]
})
names(DE_list_complete)<- contrasts$comparison
saveRDS(DE_list_complete, file.path(analysis.dir, "de_list_complete.rds"))
