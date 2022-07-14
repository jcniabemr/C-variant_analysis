#RNA Sequencing analysis of the C-300005 Fusarium venenatum A3/5 experiment. 
#Experiment consists of 154 samples with two time points 

#Create conda environment for RNA-Seq analysis 
conda create --name RNASeq
conda activate RNAseq

#Install R version 4.1.3
conda install -c conda-forge r-base=4.1.3

#Move to R
R

#Install packages 


#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")


#BiocManager::install("DESeq2")
#BiocManager::install("BiocParallel", force = TRUE)
#install.packages("ggplot2")
#BiocManager::install("Biostrings", force = TRUE)
#conda install -c conda-forge r-devtools
#install.packages("devtools")
#install.packages("data.table")
#install.packages("dplyr")
#install.packages("naturalsort")
#install.packages("tibble")
#BiocManager::install("tximport")
#install.packages("rjson")
#install.packages("readr")
#install.packages("pheatmap")
#install.packages("data.table")
#install.packages("RColorBrewer")
#install.packages("gplots")
#install.packages("ggrepel")
#conda install -c conda-forge r-wesanderson



#Set working directory 
setwd("/home/jconnell/projects/niab/johnc/RNA_seq_data")

#===============================================================================
# Load libraries
#===============================================================================


library(DESeq2)
library(BiocParallel)
register(MulticoreParam(12))
library(ggplot2)
library(Biostrings)
library(devtools) 
library(data.table)
library(dplyr)
library(naturalsort)
library(tibble)
library(tximport)
library(rjson)
library(readr)
library(pheatmap)
library(data.table)
library(RColorBrewer)
library(gplots)
library(ggrepel)
library(wesanderson)


#===============================================================================
# DeSeq analysis
#===============================================================================

# import transcript to gene mapping info
tx2gene <- read.table("DeSeq2/trans2gene.txt",header=T,sep="\t")

# import quantification files

txi.reps <- tximport(paste(list.dirs("DeSeq2", full.names=T,recursive=F),"/quant.sf",sep=""),type="salmon",tx2gene=tx2gene,txOut=T)


# get the sample names from folders

mysamples <- list.dirs("DeSeq2",full.names=F,recursive=F)

# summarise to gene level. This can be done in the tximport step (txOut=F), but is easier to understand in two steps.
txi.genes <- summarizeToGene(txi.reps,tx2gene)

names(txi.genes)


# set the sample names for txi.genes
invisible(sapply(seq(1,3), function(i) {colnames(txi.genes[[i]])<<-mysamples}))

# Get TPM tables
write.table(txi.genes,"DeSeq2/txigenes.txt",sep="\t",na="",quote=F)
write.table(txi.reps,"DeSeq2/txireps.txt",sep="\t",na="",quote=F)


# Read sample metadata
# Data is unordered as it is read in. This means data must be set into the same
# order as the samples were read into mysamples before integrating metadata and
# and read counts

unorderedColData <- read.table("DeSeq2/RNAseq_design_2.txt",header=T,sep="\t")
colData <- data.frame(unorderedColData[ order(unorderedColData$Sample),])


# Group column
#colData$Group <- paste0(colData$Condition=="Experimental", colData$Condition=="Control")
#colData$Group <- paste0(colData$Condition, colData$Timepoint)

colData$Group <- paste0(colData$Condition, colData$ID)


design <- ~ Group


dds <- DESeqDataSetFromTximport(txi.genes,colData,design)
keep <- rowSums(counts(dds)) >= 50
dds <- dds[keep,]
# Library normalisation
dds <- estimateSizeFactors(dds)
# Deseq
dds<-DESeq(dds)
resultsNames(dds)


#[1] "Intercept"                         "Group_Control.2_vs_Control.1"
#[3] "Group_Experimental.1_vs_Control.1" "Group_Experimental.2_vs_Control.1"
#[4] "Group_Experimental.2._vs_Control.2" "Group_Experimental.1_vs_Control.2"

res <- results(dds)
res
summary(res)


# MA-plot
 pdf("DeSeq2/plotMA_vst.pdf")
 plotMA(res, ylim=c(-2,2))
 dev.off()

#Adjusted p value cut off 
alpha <- 0.05 

#res= results(dds, alpha=alpha,contrast=c("Group","1","39"))
#sig.res <- subset(res,padj<=alpha)
#sig.res <- sig.res[order(sig.res$padj),]
#sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
#sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
#summary(sig.res)
#write.table(sig.res,"DeSeq2/independent_deg/TP1/1/C1_vs_WT_TP1_RES.txt",sep="\t",na="",quote=F)
#write.table(sig.res.upregulated,"DeSeq2/independent_deg/TP1/1/C1_VS_WT_TP1_UP.txt",sep="\t",na="",quote=F)
#write.table(sig.res.downregulated,"DeSeq2/independent_deg/TP1/1/C1_VS_WT_TP1_DN.txt",sep="\t",na="",quote=F)

#res= results(dds, alpha=alpha,contrast=c("Group","3","39"))
#sig.res <- subset(res,padj<=alpha)
#sig.res <- sig.res[order(sig.res$padj),]
#sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
#sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
#summary(sig.res)
#write.table(sig.res,"DeSeq2/independent_deg/TP1/2/C2_vs_WT_TP1_RES.txt",sep="\t",na="",quote=F)
#write.table(sig.res.upregulated,"DeSeq2/independent_deg/TP1/2/C2_VS_WT_TP1_UP.txt",sep="\t",na="",quote=F)
#write.table(sig.res.downregulated,"DeSeq2/independent_deg/TP1/2/C2_VS_WT_TP1_DN.txt",sep="\t",na="",quote=F)

#res= results(dds, alpha=alpha,contrast=c("Group","5","39"))
#sig.res <- subset(res,padj<=alpha)
#sig.res <- sig.res[order(sig.res$padj),]
#sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
#sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
#summary(sig.res)
#write.table(sig.res,"DeSeq2/independent_deg/TP1/3/C3_vs_WT_TP1_RES.txt",sep="\t",na="",quote=F)
#write.table(sig.res.upregulated,"DeSeq2/independent_deg/TP1/3/C3_VS_WT_TP1_UP.txt",sep="\t",na="",quote=F)
#write.table(sig.res.downregulated,"DeSeq2/independent_deg/TP1/3/C3_VS_WT_TP1_DN.txt",sep="\t",na="",quote=F)

#res= results(dds, alpha=alpha,contrast=c("Group","7","39"))
#sig.res <- subset(res,padj<=alpha)
#sig.res <- sig.res[order(sig.res$padj),]
#sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
#sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
#summary(sig.res)
#write.table(sig.res,"DeSeq2/independent_deg/TP1/4/C4_vs_WT_TP1_RES.txt",sep="\t",na="",quote=F)
#write.table(sig.res.upregulated,"DeSeq2/independent_deg/TP1/4/C4_VS_WT_TP1_UP.txt",sep="\t",na="",quote=F)
#write.table(sig.res.downregulated,"DeSeq2/independent_deg/TP1/4/C4_VS_WT_TP1_DN.txt",sep="\t",na="",quote=F)

#res= results(dds, alpha=alpha,contrast=c("Group","9","39"))
#sig.res <- subset(res,padj<=alpha)
#sig.res <- sig.res[order(sig.res$padj),]
#sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
#sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
#summary(sig.res)
#write.table(sig.res,"DeSeq2/independent_deg/TP1/5/C5_vs_WT_TP1_RES.txt",sep="\t",na="",quote=F)
#write.table(sig.res.upregulated,"DeSeq2/independent_deg/TP1/5/C5_VS_WT_TP1_UP.txt",sep="\t",na="",quote=F)
#write.table(sig.res.downregulated,"DeSeq2/independent_deg/TP1/5/C5_VS_WT_TP1_DN.txt",sep="\t",na="",quote=F)

#res= results(dds, alpha=alpha,contrast=c("Group","11","39"))
#sig.res <- subset(res,padj<=alpha)
#sig.res <- sig.res[order(sig.res$padj),]
#sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
#sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
#summary(sig.res)
#write.table(sig.res,"DeSeq2/independent_deg/TP1/6/C6_vs_WT_TP1_RES.txt",sep="\t",na="",quote=F)
#write.table(sig.res.upregulated,"DeSeq2/independent_deg/TP1/6/C6_VS_WT_TP1_UP.txt",sep="\t",na="",quote=F)
#write.table(sig.res.downregulated,"DeSeq2/independent_deg/TP1/6/C6_VS_WT_TP1_DN.txt",sep="\t",na="",quote=F)

#res= results(dds, alpha=alpha,contrast=c("Group","13","39"))
#sig.res <- subset(res,padj<=alpha)
#sig.res <- sig.res[order(sig.res$padj),]
#sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
#sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
#summary(sig.res)
#write.table(sig.res,"DeSeq2/independent_deg/TP1/7/C7_vs_WT_TP1_RES.txt",sep="\t",na="",quote=F)
#write.table(sig.res.upregulated,"DeSeq2/independent_deg/TP1/7/C7_VS_WT_TP1_UP.txt",sep="\t",na="",quote=F)
#write.table(sig.res.downregulated,"DeSeq2/independent_deg/TP1/7/C7_VS_WT_TP1_DN.txt",sep="\t",na="",quote=F)

#res= results(dds, alpha=alpha,contrast=c("Group","15","39"))
#sig.res <- subset(res,padj<=alpha)
#sig.res <- sig.res[order(sig.res$padj),]
#sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
#sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
#summary(sig.res)
#write.table(sig.res,"DeSeq2/independent_deg/TP1/8/C8_vs_WT_TP1_RES.txt",sep="\t",na="",quote=F)
#write.table(sig.res.upregulated,"DeSeq2/independent_deg/TP1/8/C8_VS_WT_TP1_UP.txt",sep="\t",na="",quote=F)
#write.table(sig.res.downregulated,"DeSeq2/independent_deg/TP1/8/C8_VS_WT_TP1_DN.txt",sep="\t",na="",quote=F)

#res= results(dds, alpha=alpha,contrast=c("Group","17","39"))
#sig.res <- subset(res,padj<=alpha)
#sig.res <- sig.res[order(sig.res$padj),]
#sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
#sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
#summary(sig.res)
#write.table(sig.res,"DeSeq2/independent_deg/TP1/9/C9_vs_WT_TP1_RES.txt",sep="\t",na="",quote=F)
#write.table(sig.res.upregulated,"DeSeq2/independent_deg/TP1/9/C9_VS_WT_TP1_UP.txt",sep="\t",na="",quote=F)
#write.table(sig.res.downregulated,"DeSeq2/independent_deg/TP1/9/C9_VS_WT_TP1_DN.txt",sep="\t",na="",quote=F)

#res= results(dds, alpha=alpha,contrast=c("Group","19","39"))
#sig.res <- subset(res,padj<=alpha)
#sig.res <- sig.res[order(sig.res$padj),]
#sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
#sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
#summary(sig.res)
#write.table(sig.res,"DeSeq2/independent_deg/TP1/10/C10_vs_WT_TP1_RES.txt",sep="\t",na="",quote=F)
#write.table(sig.res.upregulated,"DeSeq2/independent_deg/TP1/10/C10_VS_WT_TP1_UP.txt",sep="\t",na="",quote=F)
#write.table(sig.res.downregulated,"DeSeq2/independent_deg/TP1/10/C10_VS_WT_TP1_DN.txt",sep="\t",na="",quote=F)

#res= results(dds, alpha=alpha,contrast=c("Group","29","39"))
#sig.res <- subset(res,padj<=alpha)
#sig.res <- sig.res[order(sig.res$padj),]
#sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
#sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
#summary(sig.res)
#write.table(sig.res,"DeSeq2/independent_deg/TP1/15/C15_vs_WT_TP1_RES.txt",sep="\t",na="",quote=F)
#write.table(sig.res.upregulated,"DeSeq2/independent_deg/TP1/15/C15_VS_WT_TP1_UP.txt",sep="\t",na="",quote=F)
#write.table(sig.res.downregulated,"DeSeq2/independent_deg/TP1/15/C15_VS_WT_TP1_DN.txt",sep="\t",na="",quote=F)

#res= results(dds, alpha=alpha,contrast=c("Group","31","39"))
#sig.res <- subset(res,padj<=alpha)
#sig.res <- sig.res[order(sig.res$padj),]
#sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
#sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
#summary(sig.res)
#write.table(sig.res,"DeSeq2/independent_deg/TP1/16/C16_vs_WT_TP1_RES.txt",sep="\t",na="",quote=F)
#write.table(sig.res.upregulated,"DeSeq2/independent_deg/TP1/16/C16_VS_WT_TP1_UP.txt",sep="\t",na="",quote=F)
#write.table(sig.res.downregulated,"DeSeq2/independent_deg/TP1/16/C16_VS_WT_TP1_DN.txt",sep="\t",na="",quote=F)

#res= results(dds, alpha=alpha,contrast=c("Group","33","39"))
#sig.res <- subset(res,padj<=alpha)
#sig.res <- sig.res[order(sig.res$padj),]
#sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
#sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
#summary(sig.res)
#write.table(sig.res,"DeSeq2/independent_deg/TP1/17/C17_vs_WT_TP1_RES.txt",sep="\t",na="",quote=F)
#write.table(sig.res.upregulated,"DeSeq2/independent_deg/TP1/17/C17_VS_WT_TP1_UP.txt",sep="\t",na="",quote=F)
#write.table(sig.res.downregulated,"DeSeq2/independent_deg/TP1/17/C17_VS_WT_TP1_DN.txt",sep="\t",na="",quote=F)

#res= results(dds, alpha=alpha,contrast=c("Group","35","39"))
#sig.res <- subset(res,padj<=alpha)
#sig.res <- sig.res[order(sig.res$padj),]
#sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
#sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
#summary(sig.res)
#write.table(sig.res,"DeSeq2/independent_deg/TP1/18/C18_vs_WT_TP1_RES.txt",sep="\t",na="",quote=F)
#write.table(sig.res.upregulated,"DeSeq2/independent_deg/TP1/18/C18_VS_WT_TP1_UP.txt",sep="\t",na="",quote=F)
#write.table(sig.res.downregulated,"DeSeq2/independent_deg/TP1/18/C18_VS_WT_TP1_DN.txt",sep="\t",na="",quote=F)

#res= results(dds, alpha=alpha,contrast=c("Group","37","39"))
#sig.res <- subset(res,padj<=alpha)
#sig.res <- sig.res[order(sig.res$padj),]
#sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
#sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
#summary(sig.res)
#write.table(sig.res,"DeSeq2/independent_deg/TP1/19/C19_vs_WT_TP1_RES.txt",sep="\t",na="",quote=F)
#write.table(sig.res.upregulated,"DeSeq2/independent_deg/TP1/19/C19_VS_WT_TP1_UP.txt",sep="\t",na="",quote=F)
#write.table(sig.res.downregulated,"DeSeq2/independent_deg/TP1/19/C19_VS_WT_TP1_DN.txt",sep="\t",na="",quote=F)



res= results(dds, alpha=alpha,contrast=c("Group","2","40"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
summary(sig.res)
write.table(sig.res,"DeSeq2/independent_deg/TP2/1/C1_vs_WT_TP2_RES.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"DeSeq2/independent_deg/TP2/1/C1_VS_WT_TP2_UP.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"DeSeq2/independent_deg/TP2/1/C1_VS_WT_TP2_DN.txt",sep="\t",na="",quote=F)

res= results(dds, alpha=alpha,contrast=c("Group","4","40"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
summary(sig.res)
write.table(sig.res,"DeSeq2/independent_deg/TP2/2/C2_vs_WT_TP2_RES.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"DeSeq2/independent_deg/TP2/2/C2_VS_WT_TP2_UP.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"DeSeq2/independent_deg/TP2/2/C2_VS_WT_TP2_DN.txt",sep="\t",na="",quote=F)

res= results(dds, alpha=alpha,contrast=c("Group","6","40"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
summary(sig.res)
write.table(sig.res,"DeSeq2/independent_deg/TP2/3/C3_vs_WT_TP2_RES.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"DeSeq2/independent_deg/TP2/3/C3_VS_WT_TP2_UP.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"DeSeq2/independent_deg/TP2/3/C3_VS_WT_TP2_DN.txt",sep="\t",na="",quote=F)

res= results(dds, alpha=alpha,contrast=c("Group","8","40"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
summary(sig.res)
write.table(sig.res,"DeSeq2/independent_deg/TP2/4/C4_vs_WT_TP2_RES.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"DeSeq2/independent_deg/TP2/4/C4_VS_WT_TP2_UP.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"DeSeq2/independent_deg/TP2/4/C4_VS_WT_TP2_DN.txt",sep="\t",na="",quote=F)

res= results(dds, alpha=alpha,contrast=c("Group","10","40"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
summary(sig.res)
write.table(sig.res,"DeSeq2/independent_deg/TP2/5/C5_vs_WT_TP2_RES.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"DeSeq2/independent_deg/TP2/5/C5_VS_WT_TP2_UP.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"DeSeq2/independent_deg/TP2/5/C5_VS_WT_TP2_DN.txt",sep="\t",na="",quote=F)

res= results(dds, alpha=alpha,contrast=c("Group","12","40"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
summary(sig.res)
write.table(sig.res,"DeSeq2/independent_deg/TP2/6/C6_vs_WT_TP2_RES.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"DeSeq2/independent_deg/TP2/6/C6_VS_WT_TP2_UP.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"DeSeq2/independent_deg/TP2/6/C6_VS_WT_TP2_DN.txt",sep="\t",na="",quote=F)

res= results(dds, alpha=alpha,contrast=c("Group","14","40"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
summary(sig.res)
write.table(sig.res,"DeSeq2/independent_deg/TP2/7/C7_vs_WT_TP2_RES.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"DeSeq2/independent_deg/TP2/7/C7_VS_WT_TP2_UP.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"DeSeq2/independent_deg/TP2/7/C7_VS_WT_TP2_DN.txt",sep="\t",na="",quote=F)

res= results(dds, alpha=alpha,contrast=c("Group","16","40"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
summary(sig.res)
write.table(sig.res,"DeSeq2/independent_deg/TP2/8/C8_vs_WT_TP2_RES.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"DeSeq2/independent_deg/TP2/8/C8_VS_WT_TP2_UP.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"DeSeq2/independent_deg/TP2/8/C8_VS_WT_TP2_DN.txt",sep="\t",na="",quote=F)

res= results(dds, alpha=alpha,contrast=c("Group","18","40"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
summary(sig.res)
write.table(sig.res,"DeSeq2/independent_deg/TP2/9/C9_vs_WT_TP2_RES.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"DeSeq2/independent_deg/TP2/9/C9_VS_WT_TP2_UP.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"DeSeq2/independent_deg/TP2/9/C9_VS_WT_TP2_DN.txt",sep="\t",na="",quote=F)

res= results(dds, alpha=alpha,contrast=c("Group","20","40"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
summary(sig.res)
write.table(sig.res,"DeSeq2/independent_deg/TP2/10/C10_vs_WT_TP2_RES.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"DeSeq2/independent_deg/TP2/10/C10_VS_WT_TP2_UP.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"DeSeq2/independent_deg/TP2/10/C10_VS_WT_TP2_DN.txt",sep="\t",na="",quote=F)

res= results(dds, alpha=alpha,contrast=c("Group","22","40"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
summary(sig.res)
write.table(sig.res,"DeSeq2/independent_deg/TP2/11/C11_vs_WT_TP1_RES.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"DeSeq2/independent_deg/TP2/11/C11_VS_WT_TP1_UP.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"DeSeq2/independent_deg/TP2/11/C11_VS_WT_TP1_DN.txt",sep="\t",na="",quote=F)

res= results(dds, alpha=alpha,contrast=c("Group","24","40"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
summary(sig.res)
write.table(sig.res,"DeSeq2/independent_deg/TP2/12/C12_vs_WT_TP1_RES.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"DeSeq2/independent_deg/TP2/12/C12_VS_WT_TP1_UP.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"DeSeq2/independent_deg/TP2/12/C12_VS_WT_TP1_DN.txt",sep="\t",na="",quote=F)

res= results(dds, alpha=alpha,contrast=c("Group","26","40"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
summary(sig.res)
write.table(sig.res,"DeSeq2/independent_deg/TP2/13/C13_vs_WT_TP1_RES.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"DeSeq2/independent_deg/TP2/13/C13_VS_WT_TP1_UP.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"DeSeq2/independent_deg/TP2/13/C13_VS_WT_TP1_DN.txt",sep="\t",na="",quote=F)

res= results(dds, alpha=alpha,contrast=c("Group","28","40"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
summary(sig.res)
write.table(sig.res,"DeSeq2/independent_deg/TP2/14/C14_vs_WT_TP2_RES.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"DeSeq2/independent_deg/TP2/14/C14_VS_WT_TP2_UP.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"DeSeq2/independent_deg/TP2/14/C14_VS_WT_TP2_DN.txt",sep="\t",na="",quote=F)

res= results(dds, alpha=alpha,contrast=c("Group","30","40"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
summary(sig.res)
write.table(sig.res,"DeSeq2/independent_deg/TP2/15/C15_vs_WT_TP2_RES.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"DeSeq2/independent_deg/TP2/15/C15_VS_WT_TP2_UP.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"DeSeq2/independent_deg/TP2/15/C15_VS_WT_TP2_DN.txt",sep="\t",na="",quote=F)

res= results(dds, alpha=alpha,contrast=c("Group","32","40"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
summary(sig.res)
write.table(sig.res,"DeSeq2/independent_deg/TP2/16/C16_vs_WT_TP2_RES.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"DeSeq2/independent_deg/TP2/16/C16_VS_WT_TP2_UP.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"DeSeq2/independent_deg/TP2/16/C16_VS_WT_TP2_DN.txt",sep="\t",na="",quote=F)

res= results(dds, alpha=alpha,contrast=c("Group","34","40"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
summary(sig.res)
write.table(sig.res,"DeSeq2/independent_deg/TP2/17/C17_vs_WT_TP2_RES.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"DeSeq2/independent_deg/TP2/17/C17_VS_WT_TP2_UP.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"DeSeq2/independent_deg/TP2/17/C17_VS_WT_TP2_DN.txt",sep="\t",na="",quote=F)

res= results(dds, alpha=alpha,contrast=c("Group","36","40"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
summary(sig.res)
write.table(sig.res,"DeSeq2/independent_deg/TP2/18/C18_vs_WT_TP2_RES.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"DeSeq2/independent_deg/TP2/18/C18_VS_WT_TP2_UP.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"DeSeq2/independent_deg/TP2/18/C18_VS_WT_TP2_DN.txt",sep="\t",na="",quote=F)

res= results(dds, alpha=alpha,contrast=c("Group","38","40"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
summary(sig.res)
write.table(sig.res,"DeSeq2/independent_deg/TP2/19/C19_vs_WT_TP2_RES.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"DeSeq2/independent_deg/TP2/19/C19_VS_WT_TP2_UP.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"DeSeq2/independent_deg/TP2/19/C19_VS_WT_TP2_DN.txt",sep="\t",na="",quote=F)



vst1<-varianceStabilizingTransformation(dds,blind=TRUE)
write.csv(assay(vst1), file="alignment/salmon/iUK/F.venenatum/WT_minion/DeSeq2/Tables_GeneID/iUK_data_vst_T.csv")
vst2<-varianceStabilizingTransformation(dds,blind=FALSE)
write.csv(assay(vst2), file="alignment/salmon/iUK/F.venenatum/WT_minion/DeSeq2/Tables_GeneID/iUK_data_vst_F.csv")

pdf("alignment/salmon/iUK/F.venenatum/WT_minion/DeSeq2/Tables_GeneID/heatmap_vst1.pdf", width=12,height=12)
sampleDists<-dist(t(assay(vst1)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vst1$Media)
colnames(sampleDistMatrix) <- paste(vst1$Media)
colours <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap( sampleDistMatrix, trace="none", col=colours, margins=c(12,12),srtCol=45)
dev.off()

pdf("alignment/salmon/iUK/F.venenatum/WT_minion/DeSeq2/Tables_GeneID/heatmap_vst2.pdf", width=12,height=12)
sampleDists<-dist(t(assay(vst2)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vst2$Media)
colnames(sampleDistMatrix) <- paste(vst2$Media)
colours <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap( sampleDistMatrix, trace="none", col=colours, margins=c(12,12),srtCol=45)
dev.off()




































































































































































































res= results(dds, alpha=alpha,contrast=c("Group","Experimental 2","Control 2"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
summary(sig.res)
write.table(sig.res,"DeSeq2/res_experimental2_vs_control2.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"DeSeq2/upregulated_experimental2_vs_control2.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"DeSeq2/downreuglated_experimental2_vs_control2.txt",sep="\t",na="",quote=F)

res= results(dds, alpha=alpha,contrast=c("Group","Experimental 1","Experimental 2"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
summary(sig.res)
write.table(sig.res,"DeSeq2/res_experimental1_vs_experimental2.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"DeSeq2/upregulated_experimental1_vs_experimental2.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"DeSeq2/downreuglated_experimental1_vs_experimental2.txt",sep="\t",na="",quote=F)


vst1<-varianceStabilizingTransformation(dds,blind=FALSE)
write.csv(assay(vst1), file="DeSeq2/RNA_seq_data_F.csv")


vst2<-varianceStabilizingTransformation(dds,blind=TRUE)
write.csv(assay(vst2), file="DeSeq2/RNA_seq_data_T.csv")

# PCA plots 

data <- plotPCA(vst1, intgroup=c("Group"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
pca_plot<- ggplot(data, aes(PC1, PC2, color=Group)) +
geom_point(size=2) +
xlab(paste0("PC1: ",percentVar[1],"% variance")) +
ylab(paste0("PC2: ",percentVar[2],"% variance")) + geom_text_repel(aes(label=colnames(vst1))) + theme_light()
coord_fixed()
ggsave("DeSeq2/PCA_vst_false.jpeg", pca_plot, dpi=300, height=10, width=12)

data <- plotPCA(vst2, intgroup=c("Group"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
pca_plot<- ggplot(data, aes(PC1, PC2, color=Group)) +
geom_point(size=3) +
xlab(paste0("PC1: ",percentVar[1],"% variance")) +
ylab(paste0("PC2: ",percentVar[2],"% variance")) + geom_text_repel(aes(label=colnames(vst2))) + theme_light()
coord_fixed()
ggsave("DeSeq2/PCA_vst_true.jpeg", pca_plot, dpi=300, height=10, width=12)

#Cut Gene and LFC from contast files for plotting

awk '{print $1,$3}' downreuglated_experimental2_vs_control2.txt > contrast_e2_vs_c2.txt




pdf("DeSeq2/heatmap_vst1.pdf", width=12,height=12)
sampleDists<-dist(t(assay(vst1)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vst1$Condition)
colnames(sampleDistMatrix) <- paste(vst1$Timepoint)
colours <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap( sampleDistMatrix, trace="none", col=colours, margins=c(12,12),srtCol=45)
dev.off()



