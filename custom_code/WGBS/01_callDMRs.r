#Directories
data.dir <- "/home/heyj/c010-datasets/Internal/COPD/NaturComm/tWGBS/data/"

#Libraries
library(DSS)
library(bsseq)
library(doParallel)
library(ChIPseeker)
library(foreach)
library(doMC)
library(rtracklayer)

#read data
AM.BS.fit <- readRDS(file.path(data.dir,"BS.fit.fil.AM.rds"))

#run dmr analysis
pheno <- pData(AM.BS.fit )
group1AM = rownames(pData(AM.BS.fit)[pData(AM.BS.fit)$Genotype == "Tg",])
group2AM=rownames(pData(AM.BS.fit)[pData(AM.BS.fit)$Genotype == "Wt",])

#subset by chromosome
AM.BS.fit.chr<- list(NULL)
for (chr in 1:length(seqnames)) {
AM.BS.fit.chr<-  chrSelectBSseq(AM.BS.fit, seqnames =seqnames[chr], order = TRUE)
names(AM.BS.fit.chr)<- seqnames[chr]
}

#run dml test
dml.AM<- as.list()
dml.AM <-for(i in 1:length(seqnames) {
dml.AM[[i]] <- DMLtest(AM.BS.fit.chr[[i]], group1 = group1AM, group2 = group2AM, smoothing = TRUE)
}
saveRDS(dml.AM, file.path(dataAM.dir, "dml.AM.rds"))

#Combine test 
com_dml.AM<- do.call("rbind",dml.AM)
#and adjust pvalues
com_dml.AM$fdr<- p.adjust(com_dml.AM$pval, method="BH")
class(com_dml.AM)[2]<-  "DMLtest"


#Call Sig DML and DMR
sig_dmls_AM <- callDML(com_dml.AM ,delta=0.1, p.threshold=0.05)
sig_dmrs_AM <- callDMR(com_dml.AM,delta=0.1, p.threshold=0.05,
             minlen=50, minCG=3, dis.merge=50, pct.sig=0.5)
sig_dmrs_AM$direction <- ifelse(sig_dmrs_AM $diff>0, "hyper","hypo")
table(sig_dmrs_AM$direction)
saveRDS(sig_dmls_AM, file.path(data.dir , "sig_dmls_AM.rds"))
sig_dmrs_AM<- readRDS(file.path(dataAM.dir , "sig_dmrs_AM.rds"))


#Annotate DMRs
#Load TxDb file
require(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

#... and convert to granges
dmrs_gr_AM <- GRanges(
  seqnames = sig_dmrs_AM$chr,
  ranges = IRanges(start = sig_dmrs_AM$start,
                   end = sig_dmrs_AM$end
  )
)
mcols(dmrs_gr_AM ) <-sig_dmrs_AM[,4:9]

#Annotation:chipseeker
dmrs_gr_anno_AM <- annotatePeak(peak= dmrs_gr_AM, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Mm.eg.db")
dmrs_gr_anno_AM.df <-as.data.frame(dmrs_gr_anno_AM)
dmrs_gr_AM <- makeGRangesFromDataFrame(dmrs_gr_anno_AM.df, keep.extra.columns=TRUE)
saveRDS(dmrs_gr_AM, file.path(data.dir, "AM_Sig_DMR.rds"))