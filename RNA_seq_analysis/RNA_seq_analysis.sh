# RNA-Seq analysis performed on Fusarium venenatum C-variant strains (19) and WT strain. 


# Files were copied from NIAB HPC to Crop Diversity HPC. 
scp -r /archives/2022_camb_general/X204SC22021418-Z01-F001/raw_data/ jconnell@gruffalo.cropdiversity.ac.uk:/home/jconnell/projects/niab/johnc/RNA_seq_data
scp -r /home/connellj/git_repos/* jconnell@gruffalo.cropdiversity.ac.uk:/home/jconnell/git_repos
scp -r /home/connellj/miniconda2/bin/fastq-mcf jconnell@gruffalo.cropdiversity.ac.uk:/home/jconnell/miniconda3/bin



# For simplicity files were renamed to contain only the strain/timepoint/repeat 

for file in projects/niab/johnc/rna_seq_data/raw_data/*; do
  mv "$file" projects/niab/johnc/rna_seq_data/raw_data/"${file#*R}"
done

# For isolates with 3 didgets a prefix "0" assed to make 3rd and 4th number timepoint and repeat  
 
for file in 122; do 
  mv "$file" "0${file}" 
done


# Perform qc on RNA-Seq data with FastQC

for strain in 0111 0112 0113 0114 0121 0122 0123 0124 0211 0212 0213 0214 0221 0222 0223 0224 0311 0312 0313 0314 0321 0322 0323 0324 0411 0412 0413 0414 0421 0422 0423 0424 0511 0512 0513 0514 0521 0522 0523 0611 0612 0613 0614 0621 0622 0623 0624 0711 0712 0713 0714 0722 0723 0724 0811 0812 0813 0814 0822 0823 0824 0911 0912 0913 0914 0921 0922 0923 0924 1011 1012 1013 1014 1021 1023 1024 1112 1113 1114 1122 1123 1124 1211 1212 1214 1221 1222 1223 1224 1311 1312 1313 1314 1321 1322 1323 1324 1411 1412 1413 1414 1422 1423 1424 1511 1512 1513 1514 1521 1522 1523 1524 1611 1612 1613 1614 1622 1623 1624 1711 1712 1713 1714 1721 1722 1723 1724 1811 1812 1814 1821 1822 1823 1824 1911 1912 1913 1914 1921 1923 1924 WT11 WT12 WT13 WT14 WT21 WT22 WT23 WT24; do
  for RawData in $(ls projects/niab/johnc/RNA_seq_data/raw_data/$strain/*.fq.gz); do
     echo $RawData  
     Repeat=$(echo $RawData | rev | cut -d '/' -f2 | sed -r 's/.{3}$//g')
     echo $Repeat
     Timepoint=$(echo $RawData | rev | cut -d '/' -f2 | rev | sed -r 's/.{1}$//g' | sed -r 's/.{2}//g')
     echo $Timepoint
     OutDir=projects/niab/johnc/RNA_seq_data/pre_QC/$strain/
     mkdir -p $OutDir
     ProgDir=git_repos/emr_repos/mycoprotein_mutant_analysis/genome_assembly
     sbatch -p long $ProgDir/fastqc.sh $RawData $OutDir
  done
done  


   
# Run fastq-mcf

for strain in 0111 0112 0113 0114 0121 0122 0123 0124 0211 0212 0213 0214 0221 0222 0223 0224 0311 0312 0313 0314 0321 0322 0323 0324 0411 0412 0413 0414 0421 0422 0423 0424 0511 0512 0513 0514 0521 0522 0523 0611 0612 0613 0614 0621 0622 0623 0624 0711 0712 0713 0714 0722 0723 0724 0811 0812 0813 0814 0822 0823 0824 0911 0912 0913 0914 0921 0922 0923 0924 1011 1012 1013 1014 1021 1023 1024 1112 1113 1114 1122 1123 1124 1211 1212 1214 1221 1222 1223 1224 1311 1312 1313 1314 1321 1322 1323 1324 1411 1412 1413 1414 1422 1423 1424 1511 1512 1513 1514 1521 1522 1523 1524 1611 1612 1613 1614 1622 1623 1624 1711 1712 1713 1714 1721 1722 1723 1724 1811 1812 1814 1821 1822 1823 1824 1911 1912 1913 1914 1921 1923 1924 WT11 WT12 WT13 WT14 WT21 WT22 WT23 WT24; do
  for file in $(ls -d projects/niab/johnc/RNA_seq_data/raw_data/$strain); do
      FileF=$(ls $file/*_1.fq.gz)
      FileR=$(ls $file/*_2.fq.gz)
      echo $FileF
      echo $FileR
      OutDir=projects/niab/johnc/RNA_seq_data/trimmed/$strain/
      mkdir -p $OutDir
      IluminaAdapters=/home/jconnell/projects/niab/johnc/RNA_seq_data/illumina_full_adapters.fa
      ProgDir=git_repos/emr_repos/Fv_C-variants/RNA_seq_analysis
      sbatch -p long $ProgDir/fastq-mcf_himem.sh $FileF $FileR $IluminaAdapters DNA $OutDir    
  done
done 


# Perform qc on RNA-Seq data with FastQC


for strain in 0111 0112 0113 0114 0121 0122 0123 0124 0211 0212 0213 0214 0221 0222 0223 0224 0311 0312 0313 0314 0321 0322 0323 0324 0411 0412 0413 0414 0421 0422 0423 0424 0511 0512 0513 0514 0521 0522 0523 0611 0612 0613 0614 0621 0622 0623 0624 0711 0712 0713 0714 0722 0723 0724 0811 0812 0813 0814 0822 0823 0824 0911 0912 0913 0914 0921 0922 0923 0924 1011 1012 1013 1014 1021 1023 1024 1112 1113 1114 1122 1123 1124 1211 1212 1214 1221 1222 1223 1224 1311 1312 1313 1314 1321 1322 1323 1324 1411 1412 1413 1414 1422 1423 1424 1511 1512 1513 1514 1521 1522 1523 1524 1611 1612 1613 1614 1622 1623 1624 1711 1712 1713 1714 1721 1722 1723 1724 1811 1812 1814 1821 1822 1823 1824 1911 1912 1913 1914 1921 1923 1924 WT11 WT12 WT13 WT14 WT21 WT22 WT23 WT24; do
  for RawData in $(ls projects/niab/johnc/RNA_seq_data/trimmed/$strain/*/*.fq.gz); do
     echo $RawData  
     Repeat=$(echo $RawData | rev | cut -d '/' -f2 | sed -r 's/.{3}$//g')
     echo $Repeat
     Timepoint=$(echo $RawData | rev | cut -d '/' -f2 | rev | sed -r 's/.{1}$//g' | sed -r 's/.{2}//g')
     echo $Timepoint
     OutDir=projects/niab/johnc/RNA_seq_data/post_QC/$strain/
     mkdir -p $OutDir
     ProgDir=git_repos/emr_repos/mycoprotein_mutant_analysis/genome_assembly
     sbatch -p long $ProgDir/fastqc.sh $RawData $OutDir
  done
done  


## Decontamination of rRNA reads in RNAseq data

for Strain in 0111 0112 0113 0114 0121 0122 0123 0124 0211 0212 0213 0214 0221 0222 0223 0224 0311 0312 0313 0314 0321 0322 0323 0324 0411 0412 0413 0414 0421 0422 0423 0424 0511 0512 0513 0514 0521 0522 0523 0611 0612 0613 0614 0621 0622 0623 0624 0711 0712 0713 0714 0722 0723 0724 0811 0812 0813 0814 0822 0823 0824 0911 0912 0913 0914 0921 0922 0923 0924 1011 1012 1013 1014 1021 1023 1024 1112 1113 1114 1122 1123 1124 1211 1212 1214 1221 1222 1223 1224 1311 1312 1313 1314 1321 1322 1323 1324 1411 1412 1413 1414 1422 1423 1424 1511 1512 1513 1514 1521 1522 1523 1524 1611 1612 1613 1614 1622 1623 1624 1711 1712 1713 1714 1721 1722 1723 1724 1811 1812 1814 1821 1822 1823 1824 1911 1912 1913 1914 1921 1923 1924 WT11 WT12 WT13 WT14 WT21 WT22 WT23 WT24; do
    for RawData in /home/jconnell/projects/niab/johnc/RNA_seq_data/trimmed/$Strain; do
        FileF=$RawData/F/*.fq.gz
        FileR=$RawData/R/*.fq.gz
        echo $FileF
        echo $FileR
        Ref=/home/jconnell/projects/niab/johnc/RNA_seq_data/ribokmers.fa.gz
        outdir=/home/jconnell/projects/niab/johnc/RNA_seq_data/decontaminated/"$Strain"
        mkdir -p $outdir
        ProgDir=/home/jconnell/git_repos/emr_repos/Fv_C-variants/RNA_seq_analysis  
        sbatch -p long $ProgDir/bbduk.sh $Ref $FileF $FileR $Strain $outdir
    done
done



## Salmon 

```bash
conda activate salmon
cd /projects/fusarium_venenatum


# Samples name corrected from Novogene form

Transcriptome=/home/jconnell/projects/niab/johnc/RNA_seq_data/final_genes_appended_renamed.cdna.fasta
for Strain in 0111 0112 0113 0114 0121 0122 0123 0124 0211 0212 0213 0214 0221 0222 0223 0224 0311 0312 0313 0314 0321 0322 0323 0324 0411 0412 0413 0414 0421 0422 0423 0424 0511 0512 0513 0514 0521 0522 0523 0611 0612 0613 0614 0621 0622 0623 0624 0711 0712 0713 0714 0722 0723 0724 0811 0812 0813 0814 0822 0823 0824 0911 0912 0913 0914 0921 0922 0923 0924 1011 1012 1013 1014 1021 1023 1024 1112 1113 1114 1122 1123 1124 1211 1212 1214 1221 1222 1223 1224 1311 1312 1313 1314 1321 1322 1323 1324 1411 1412 1413 1414 1422 1423 1424 1511 1512 1513 1514 1521 1522 1523 1524 1611 1612 1613 1614 1622 1623 1624 1711 1712 1713 1714 1721 1722 1723 1724 1811 1812 1814 1821 1822 1823 1824 1911 1912 1913 1914 1921 1923 1924 WT11 WT12 WT13 WT14 WT21 WT22 WT23 WT24; do
    for RawData in /home/jconnell/projects/niab/johnc/RNA_seq_data/decontaminated/$Strain; do
        FileF=$RawData/F/*cleaned.fq.gz
        FileR=$RawData/R/*cleaned.fq.gz
        echo $FileF
        echo $FileR
        OutDir=/home/jconnell/projects/niab/johnc/RNA_seq_data/salmon/$Strain
        mkdir -p $OutDir
        ProgDir=/home/jconnell/git_repos/emr_repos/Fv_C-variants/RNA_seq_analysis
        sbatch -p long $ProgDir/salmon.sh $Transcriptome $FileF $FileR $OutDir
    done
done



#Convert Salmon quasi-quanitifcations to gene counts using an awk script:


for Strain in 0111 0112 0113 0114 0121 0122 0123 0124 0211 0212 0213 0214 0221 0222 0223 0224 0311 0312 0313 0314 0321 0322 0323 0324 0411 0412 0413 0414 0421 0422 0423 0424 0511 0512 0513 0514 0521 0522 0523 0611 0612 0613 0614 0621 0622 0623 0624 0711 0712 0713 0714 0722 0723 0724 0811 0812 0813 0814 0822 0823 0824 0911 0912 0913 0914 0921 0922 0923 0924 1011 1012 1013 1014 1021 1023 1024 1112 1113 1114 1122 1123 1124 1211 1212 1214 1221 1222 1223 1224 1311 1312 1313 1314 1321 1322 1323 1324 1411 1412 1413 1414 1422 1423 1424 1511 1512 1513 1514 1521 1522 1523 1524 1611 1612 1613 1614 1622 1623 1624 1711 1712 1713 1714 1721 1722 1723 1724 1811 1812 1814 1821 1822 1823 1824 1911 1912 1913 1914 1921 1923 1924 WT11 WT12 WT13 WT14 WT21 WT22 WT23 WT24; do
  for File in $(ls /home/jconnell/projects/niab/johnc/RNA_seq_data/salmon/$Strain/quant.sf | head -n1); do
      outdir=/home/jconnell/projects/niab/johnc/RNA_seq_data/gene_count/$Strain
      mkdir -p $outdir
      cat $File | awk -F"\t" '{c=$1;sub(".t.*","",$1);print c,$1}' OFS="\t" > trans2gene.txt
      mv trans2gene.txt $outdir
  done
done 



mkdir -p alignment/salmon/Fvenenatum_CarbonRNAseq/F.venenatum/WT_minion/DeSeq2
#This command creates a two column file with transcript_id and gene_id.
for File in $(ls alignment/salmon/Fvenenatum_CarbonRNAseq/F.venenatum/WT_minion/*/*/*/quant.sf | head -n1); do
cat $File | awk -F"\t" '{c=$1;sub(".t.*","",$1);print c,$1}' OFS="\t" > alignment/salmon/Fvenenatum_CarbonRNAseq/F.venenatum/WT_minion/DeSeq2/trans2gene.txt
done

#This command creates a two column file with transcript_id.
#for File in $(ls alignment/salmon/*/Hg199/*/*/*/quant.sf | head -n1); do
#cat $File | awk -F"\t" '{c=$1;sub("\*","",$1);print c,$1}' OFS="\t" > alignment/salmon/N.ditissima/Hg199/DeSeq2/trans2gene3.txt
#done


for Data in 0111 0112 0113 0114 0121 0122 0123 0124 0211 0212 0213 0214 0221 0222 0223 0224 0311 0312 0313 0314 0321 0322 0323 0324 0411 0412 0413 0414 0421 0422 0423 0424 0511 0512 0513 0514 0521 0522 0523 0611 0612 0613 0614 0621 0622 0623 0624 0711 0712 0713 0714 0722 0723 0724 0811 0812 0813 0814 0822 0823 0824 0911 0912 0913 0914 0921 0922 0923 0924 1011 1012 1013 1014 1021 1023 1024 1112 1113 1114 1122 1123 1124 1211 1212 1214 1221 1222 1223 1224 1311 1312 1313 1314 1321 1322 1323 1324 1411 1412 1413 1414 1422 1423 1424 1511 1512 1513 1514 1521 1522 1523 1524 1611 1612 1613 1614 1622 1623 1624 1711 1712 1713 1714 1721 1722 1723 1724 1811 1812 1814 1821 1822 1823 1824 1911 1912 1913 1914 1921 1923 1924 WT11 WT12 WT13 WT14 WT21 WT22 WT23 WT24 1421 0812 1121 0524 1022; do  
  for File in projects/niab/johnc/RNA_seq_data/salmon/$Data/quant.sf; do
    location=projects/niab/johnc/RNA_seq_data/quantification/$Data
    mkdir -p $location
    cp $File $location
  done   
done 



# Put files in a convenient location for DeSeq.
# Analysis was not performed on Apple control samples.

for File in $(ls alignment/salmon/Fvenenatum_CarbonRNAseq/F.venenatum/WT_minion/*/*/*/quant.sf); do
Prefix=$(echo $File | cut -f8 -d '/' --output-delimiter '_')
  mkdir -p alignment/salmon/Fvenenatum_CarbonRNAseq/F.venenatum/WT_minion/DeSeq2/$Prefix
  cp $PWD/$File alignment/salmon/Fvenenatum_CarbonRNAseq/F.venenatum/WT_minion/DeSeq2/$Prefix/quant.sf
  # rm alignment/salmon/DeSeq2/$Prefix/quant.sf
done
```


#for strain in 0111 0112 0113 0114 0121 0122 0123 0124 0211 0212 0213 0214 0221 0222 0223 0224 0311 0312 0313 0314 0321 0322 0323 0324 0411 0412 0413 0414 0421 0422 0423 0424 0511 0512 0513 0514 0521 0522 0523 0611 0612 0613 0614 0621 0622 0623 0624 0711 0712 0713 0714 0722 0723 0724 0811 0812 0813 0814 0822 0823 0824 0911 0912 0913 0914 0921 0922 0923 0924 1011 1012 1013 1014 1021 1023 1024 1112 1113 1114 1122 1123 1124 1211 1212 1214 1221 1222 1223 1224 1311 1312 1313 1314 1321 1322 1323 1324 1411 1412 1413 1414 1422 1423 1424 1511 1512 1513 1514 1521 1522 1523 1524 1611 1612 1613 1614 1622 1623 1624 1711 1712 1713 1714 1721 1722 1723 1724 1811 1812 1814 1821 1822 1823 1824 1911 1912 1913 1914 1921 1923 1924 WT11 WT12 WT13 WT14 WT21 WT22 WT23 WT24 1421 0812 1121 0524 1022; do 
 # for File in projects/niab/johnc/RNA_seq_data/salmon/$strain/quant.sf; do
 #   new_dir=projects/niab/johnc/RNA_seq_data/salmon_counts
 ##   mkdir -p $new_dir
  #  cp $File $new_dir/"$strain"_quant.sf 
  #  done 
#done     

# Differential expression with DeSeq


```bash
/projects/software/R-3.6.1/bin/R
```

```R
setwd("/projects/fusarium_venenatum")

#===============================================================================
#       Load libraries
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


# Load data from SALMON quasi mapping

# Analysis in DeSeq2 folder include all samples. Deseq_v2 does not include C0T0 samples.

# import transcript to gene mapping info
tx2gene <- read.table("alignment/salmon/Fvenenatum_CarbonRNAseq/F.venenatum/WT_minion/DeSeq2/trans2gene.txt",header=T,sep="\t")

# import quantification files

txi.reps <- tximport(paste(list.dirs("/home/jconnell/projects/niab/johnc/RNA_seq_data/DeSeq2", full.names=T,recursive=F),"/quant.sf",sep=""),type="salmon",tx2gene=tx2gene,txOut=T)
# No C0T0
#txi.reps <- tximport(paste(list.dirs("alignment/salmon/Fvenenatum_CarbonRNAseq/F.venenatum/WT_minion/DeSeq2_v2", full.names=T,recursive=F),"/quant.sf",sep=""),type="salmon",tx2gene=tx2gene,txOut=T)

# get the sample names from the folders

mysamples <- list.dirs("alignment/salmon/Fvenenatum_CarbonRNAseq/F.venenatum/WT_minion/DeSeq2",full.names=F,recursive=F)
mysamples <- list.dirs("alignment/salmon/Fvenenatum_CarbonRNAseq/F.venenatum/WT_minion/DeSeq2_v2",full.names=F,recursive=F)

# summarise to gene level. This can be done in the tximport step (txOut=F), but is easier to understand in two steps.
txi.genes <- summarizeToGene(txi.reps,tx2gene)

names(txi.genes)

# set the sample names for txi.genes
invisible(sapply(seq(1,3), function(i) {colnames(txi.genes[[i]])<<-mysamples}))

# Get TPM tables
write.table(txi.genes,"alignment/salmon/Fvenenatum_CarbonRNAseq/F.venenatum/WT_minion/DeSeq2/txigenes.txt",sep="\t",na="",quote=F)
write.table(txi.reps,"alignment/salmon/Fvenenatum_CarbonRNAseq/F.venenatum/WT_minion/DeSeq2/txireps.txt",sep="\t",na="",quote=F)

# PCA TPMs
pca(experiment.table="alignment/salmon/Fvenenatum_CarbonRNAseq/F.venenatum/WT_minion/DeSeq2/txigenes_TPMonly.txt", type="TPM",
legend.position="topleft", covariatesInNames=FALSE, samplesName=TRUE,
principal.components=c(1,2), pdf = TRUE,
output.folder=getwd())

pca(experiment.table="alignment/salmon/Fvenenatum_CarbonRNAseq/F.venenatum/WT_minion/DeSeq2/txireps_TPMonly.txt", type="TPM",
legend.position="topleft", covariatesInNames=FALSE, samplesName=TRUE,
principal.components=c(1,2), pdf = TRUE,
output.folder=getwd())

# Read sample metadata
# Data is unordered as it is read in. This means data must be set into the same
# order as the samples were read into mysamples before integrating metadata and
# and read counts

unorderedColData <- read.table("alignment/salmon/Fvenenatum_CarbonRNAseq/F.venenatum/WT_minion/DeSeq2/FvenCarbon_RNAseq_design3.txt",header=T,sep="\t")
colData <- data.frame(unorderedColData[ order(unorderedColData$Sample),])

# Group column
colData$Group <- paste0(colData$Condition,'_', colData$Timepoint)

# Group column
colData$Group2 <- paste0(colData$Condition.1,'_', colData$Timepoint)

# Define the DESeq 'GLM' model
design <- ~ Group
dds <- DESeqDataSetFromTximport(txi.genes,colData,design)

keep <- rowSums(counts(dds)) >= 50
dds <- dds[keep,]

# Library normalisation
dds <- estimateSizeFactors(dds)

# Deseq
dds<-DESeq(dds) paralell??

resultsNames(dds)
###
 [1] "Intercept"            "Group_C1_T1_vs_C0_T0" "Group_C1_T2_vs_C0_T0"
 [4] "Group_C1_T3_vs_C0_T0" "Group_C1_T4_vs_C0_T0" "Group_C1_T5_vs_C0_T0"
 [7] "Group_C1_T6_vs_C0_T0" "Group_C1_T7_vs_C0_T0" "Group_C2_T1_vs_C0_T0"
[10] "Group_C2_T2_vs_C0_T0" "Group_C2_T3_vs_C0_T0" "Group_C2_T4_vs_C0_T0"
[13] "Group_C2_T5_vs_C0_T0" "Group_C2_T6_vs_C0_T0" "Group_C2_T7_vs_C0_T0"
[16] "Group_C3_T1_vs_C0_T0" "Group_C3_T2_vs_C0_T0" "Group_C3_T3_vs_C0_T0"
[19] "Group_C3_T4_vs_C0_T0" "Group_C3_T5_vs_C0_T0" "Group_C3_T6_vs_C0_T0"
[22] "Group_C3_T7_vs_C0_T0" "Group_C4_T1_vs_C0_T0" "Group_C4_T2_vs_C0_T0"
[25] "Group_C4_T3_vs_C0_T0" "Group_C4_T4_vs_C0_T0" "Group_C4_T5_vs_C0_T0"
[28] "Group_C4_T6_vs_C0_T0" "Group_C4_T7_vs_C0_T0"
###


# Exploring and exporting results

res <- results(dds)
res
summary(res)

alpha <- 0.05

res= results(dds, alpha=alpha,contrast=c("Condition","C1","C0"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
summary(sig.res)
###
out of 4192 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)       : 2425, 58%
LFC < 0 (down)     : 1767, 42%
outliers [1]       : 0, 0%
low counts [2]     : 0, 0%
(mean count < 0)
###
write.table(sig.res,"C1_vs_C0.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"C1_vs_C0_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"C1_vs_C0_down.txt",sep="\t",na="",quote=F)

res= results(dds, alpha=alpha,contrast=c("Condition","C2","C0"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
summary(sig.res)
###
out of 4192 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)       : 2425, 58%
LFC < 0 (down)     : 1767, 42%
outliers [1]       : 0, 0%
low counts [2]     : 0, 0%
(mean count < 0)
###
write.table(sig.res,"C2_vs_C0.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"C2_vs_C0_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"C2_vs_C0_down.txt",sep="\t",na="",quote=F)

res= results(dds, alpha=alpha,contrast=c("Condition","C3","C0"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
summary(sig.res)
###
out of 4192 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)       : 2425, 58%
LFC < 0 (down)     : 1767, 42%
outliers [1]       : 0, 0%
low counts [2]     : 0, 0%
(mean count < 0)
###
write.table(sig.res,"C3_vs_C0.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"C3_vs_C0_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"C3_vs_C0_down.txt",sep="\t",na="",quote=F)

res= results(dds, alpha=alpha,contrast=c("Condition","C4","C0"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
summary(sig.res)
###
out of 4192 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)       : 2425, 58%
LFC < 0 (down)     : 1767, 42%
outliers [1]       : 0, 0%
low counts [2]     : 0, 0%
(mean count < 0)
###
write.table(sig.res,"C4_vs_C0.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"C4_vs_C0_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"C4_vs_C0_down.txt",sep="\t",na="",quote=F)


# Sample Distances

vst1<-varianceStabilizingTransformation(dds)
vst2<-vst(dds,blind=FALSE)
vst3<-vst(dds,blind=TRUE)
pdf("alignment/salmon/Fvenenatum_CarbonRNAseq/F.venenatum/WT_minion/DeSeq2/heatmap_vst.pdf", width=12,height=12)
sampleDists<-dist(t(assay(vst)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vst$Condition)
colnames(sampleDistMatrix) <- paste(vst$Timepoint)
colours <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap( sampleDistMatrix, trace="none", col=colours, margins=c(12,12),srtCol=45)
#heatmap( sampleDistMatrix,
#  trace="none",  # turns off trace lines inside the heat map
#  col=colours, # use on color palette defined earlier
#  margins=c(12,12), # widens margins around plot
#  srtCol=45,
#  srtCol=45)
dev.off()

# Sample distances measured with rlog transformation:
 rld <- rlog(dds)
# pdf("alignment/salmon/Fvenenatum_CarbonRNAseq/F.venenatum/WT_minion/DeSeq2/heatmap_rld.pdf")
# sampleDists <- dist(t(assay(rld)))
# sampleDistMatrix <- as.matrix( sampleDists )
# rownames(sampleDistMatrix) <- paste(rld$Cultivar)
# colnames(sampleDistMatrix) <- paste(rld$Timepoint)
# colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
# heatmap( sampleDistMatrix, trace="none", col=colours, margins=c(12,12),srtCol=45)
# dev.off()

# # MA-plot
# pdf("alignment/salmon/N.ditissima/Hg199/DeSeq2/plotMA_vst.pdf")
# plotMA(res, ylim=c(-2,2))
# dev.off()

# Plot counts
# pdf("alignment/salmon/N.ditissima/Hg199/DeSeq2/plotcounts_dds.pdf")
# plotCounts(dds, gene=which.min(res$padj), intgroup="Cultivar")
# dev.off()

# pdf("alignment/salmon/N.ditissima/Hg199/DeSeq2/plotcounts2_dds.pdf")
# plotCounts(dds, gene=which.min(res$padj), intgroup=c("Cultivar","Timepoint"))
# dev.off()

# PCA plots
pdf("alignment/salmon/Fvenenatum_CarbonRNAseq/F.venenatum/WT_minion/DeSeq2/PCA_vst_group_filtered.pdf")
plotPCA(vst,intgroup=c("Group"))
dev.off()

pdf("alignment/salmon/Fvenenatum_CarbonRNAseq/F.venenatum/WT_minion/DeSeq2_v2/PCA_vst_FALSE.pdf")
plotPCA(vst2,intgroup=c("Group"))
dev.off()

pdf("alignment/salmon/Fvenenatum_CarbonRNAseq/F.venenatum/WT_minion/DeSeq2_v2/PCA_vst_TRUE.pdf")
plotPCA(vst3,intgroup=c("Group"))
dev.off()

#Plot using rlog transformation:
pdf("alignment/salmon/Fvenenatum_CarbonRNAseq/F.venenatum/WT_minion/DeSeq2/PCA_rld_group.pdf")
plotPCA(rld,intgroup=c("Timepoint"))
dev.off()

#Plot using rlog transformation, showing sample names:

data <- plotPCA(vst, intgroup=c("Condition.1"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
pca_plot<- ggplot(data, aes(PC1, PC2, color=Condition.1)) +
geom_point(size=3) +
xlab(paste0("PC1: ",percentVar[1],"% variance")) +
ylab(paste0("PC2: ",percentVar[2],"% variance")) + geom_text_repel(aes(label=colnames(vst)))
coord_fixed()
ggsave("alignment/salmon/Fvenenatum_CarbonRNAseq/F.venenatum/WT_minion/DeSeq2/PCA_sample_names4.pdf", pca_plot, dpi=300, height=10, width=12)

data <- plotPCA(rld, intgroup=c("Condition","Timepoint"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
pca_plot<- ggplot(data, aes(PC1, PC2, color=Condition)) +
geom_point(size=3) +
xlab(paste0("PC1: ",percentVar[1],"% variance")) +
ylab(paste0("PC2: ",percentVar[2],"% variance")) + geom_text_repel(aes(label=colnames(rld)))
coord_fixed()
ggsave("alignment/salmon/Fvenenatum_CarbonRNAseq/F.venenatum/WT_minion/DeSeq2/PCA_rld_sample_names.pdf", pca_plot, dpi=300, height=10, width=12)


# Gene clustering plots

topVarGenes <-head(order(rowVars(assay(rld)),decreasing=TRUE),50)
pdf("alignment/salmon/N.ditissima/Hg199/DeSeq2/heatmap_new.pdf")
heatmap.2(assay(rld)[topVarGenes,],ColSideColors=c("grey","dodgerblue")[ rld$Timepoint ],scale='row',trace="none",dendrogram="column",col=colorRampPalette(rev(brewer.pal(9,"RdBu")))(255))
dev.off()

topVarGenes <-head(order(rowVars(assay(rld)),decreasing=TRUE),50)
pdf("alignment/salmon/N.ditissima/Hg199/DeSeq2/heatmap_new2.pdf")
colors <- colorRampPalette( rev(brewer.pal(9, "PuOr")) )(255)
sidecols <- c("grey","dodgerblue")[ rld$Cultivar ]
mat <- assay(rld)[ topVarGenes, ]
#mat <- mat - rowMeans(mat)
#colnames(mat) <- paste0(rld$Cultivar,"-",rld$Timepoint)
heatmap.2(mat, trace="none", col=colors, dendrogram="column",ColSideColors=sidecols,labRow=TRUE, mar=c(10,2), scale="row")
dev.off()
```

```r
# Done but not needed
# Make a table of raw counts, normalised counts and fpkm values:
raw_counts <- data.frame(counts(dds, normalized=FALSE))
colnames(raw_counts) <- paste(colData$Sample)
write.table(raw_counts,"alignment/salmon/Fvenenatum_CarbonRNAseq/F.venenatum/WT_minion/DeSeq2/raw_counts.txt",sep="\t",na="",quote=F)
norm_counts <- data.frame(counts(dds, normalized=TRUE))
colnames(norm_counts) <- paste(colData$Sample)
write.table(norm_counts,"alignment/salmon/Fvenenatum_CarbonRNAseq/F.venenatum/WT_minion/DeSeq2/normalised_counts.txt",sep="\t",na="",quote=F)

# robust may be better set at false to normalise based on total counts rather than 'library normalisation factors'
fpkm_counts <- data.frame(fpkm(dds, robust = TRUE))
colnames(fpkm_counts) <- paste(colData$Sample)
write.table(fpkm_counts,"alignment/salmon/Fvenenatum_CarbonRNAseq/F.venenatum/WT_minion/DeSeq2/fpkm_norm_counts.txt",sep="\t",na="",quote=F)
fpkm_counts <- data.frame(fpkm(dds, robust = FALSE))
colnames(fpkm_counts) <- paste(colData$Sample)
write.table(fpkm_counts,"alignment/salmon/Fvenenatum_CarbonRNAseq/F.venenatum/WT_minion/DeSeq2/fpkm_counts.txt",sep="\t",na="",quote=F)


pca(experiment.table="raw_counts.txt", type="counts",
      legend.position="topleft", covariatesInNames=FALSE, samplesName=TRUE,
      principal.components=c(1,2), pdf = TRUE,
      output.folder=getwd())


pca(experiment.table="normalised_counts.txt", type="counts",
      legend.position="topleft", covariatesInNames=FALSE, samplesName=TRUE,
      principal.components=c(1,2), pdf = TRUE,
      output.folder=getwd())

pca(experiment.table="fpkm_counts.txt", type="FPKM",
      legend.position="topleft", covariatesInNames=FALSE, samplesName=TRUE,
      principal.components=c(1,2), pdf = TRUE,
      output.folder=getwd())

pca(experiment.table="fpkm_norm_counts.txt", type="FPKM",
      legend.position="topleft", covariatesInNames=FALSE, samplesName=TRUE,
      principal.components=c(1,2), pdf = TRUE,
      output.folder=getwd())


write.csv(vst, file="vst_all.csv")
write.csv(assay(vst), file="vst_all.csv")
```

```r
# For cosistently expressed genes

topVarGenes <- head( order( rowVars( assay(vst) ), decreasing=FALSE ), 1000 ) # decreasing=TRUE for most variable genes
mat <- assay(vst)[ topVarGenes, ]
write.table(mat,"Topvar2.txt",sep="\t",na="",quote=F)

# Combine both datasets and look shared genes
T1<-read.table("Topvar1.txt",header=T,sep="\t")
T2<-read.table("Topvar2.txt",header=T,sep="\t")
T3<-merge(T1,T2, by.x="ID",by.y="ID",all.x=TRUE,all.y=TRUE) # Print all
T3<-merge(T1,T2, by.x="ID",by.y="ID",all.x=FALSE,all.y=FALSE) # Shared only
write.table(T3,"common.txt",sep="\t",na="",quote=F)
```





# Corrected data


Convert Salmon quasi-quanitifcations to gene counts using an awk script:

```bash
mkdir -p alignment/salmon/Fvenenatum_CarbonRNAseq_CORRECTED/F.venenatum/WT_minion/corrected/DeSeq2
#This command creates a two column file with transcript_id and gene_id.
for File in $(ls alignment/salmon/Fvenenatum_CarbonRNAseq_CORRECTED/F.venenatum/WT_minion/corrected/*/*/quant.sf | head -n1); do
cat $File | awk -F"\t" '{c=$1;sub(".t.*","",$1);print c,$1}' OFS="\t" > alignment/salmon/Fvenenatum_CarbonRNAseq_CORRECTED/F.venenatum/WT_minion/corrected/DeSeq2/trans2gene.txt
done

#This command creates a two column file with transcript_id.
#for File in $(ls alignment/salmon/*/Hg199/*/*/*/quant.sf | head -n1); do
#cat $File | awk -F"\t" '{c=$1;sub("\*","",$1);print c,$1}' OFS="\t" > alignment/salmon/N.ditissima/Hg199/DeSeq2/trans2gene3.txt
#done

# Put files in a convenient location for DeSeq.
# Analysis was not performed on Apple control samples.

for File in $(ls alignment/salmon/Fvenenatum_CarbonRNAseq_CORRECTED/F.venenatum/WT_minion/corrected/*/*/quant.sf); do
Prefix=$(echo $File | cut -f8 -d '/' --output-delimiter '_')
mkdir -p alignment/salmon/Fvenenatum_CarbonRNAseq_CORRECTED/F.venenatum/WT_minion/corrected/DeSeq2/$Prefix
cp $PWD/$File alignment/salmon/Fvenenatum_CarbonRNAseq_CORRECTED/F.venenatum/WT_minion/corrected/DeSeq2/$Prefix/quant.sf
# rm alignment/salmon/DeSeq2/$Prefix/quant.sf
done
```

# Differential expression with DeSeq


```bash
/projects/software/R-3.6.1/bin/R
```

```R
setwd("/projects/fusarium_venenatum")

# Load libraries

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


# Load data from SALMON quasi mapping

# import transcript to gene mapping info
tx2gene <- read.table("alignment/salmon/Fvenenatum_CarbonRNAseq_CORRECTED/F.venenatum/WT_minion/corrected/DeSeq2/trans2gene.txt",header=T,sep="\t")

# import quantification files

txi.reps <- tximport(paste(list.dirs("alignment/salmon/Fvenenatum_CarbonRNAseq_CORRECTED/F.venenatum/WT_minion/corrected/DeSeq2", full.names=T,recursive=F),"/quant.sf",sep=""),type="salmon",tx2gene=tx2gene,txOut=T)

# get the sample names from the folders

mysamples <- list.dirs("alignment/salmon/Fvenenatum_CarbonRNAseq_CORRECTED/F.venenatum/WT_minion/corrected/DeSeq2",full.names=F,recursive=F)

# summarise to gene level. This can be done in the tximport step (txOut=F), but is easier to understand in two steps.
txi.genes <- summarizeToGene(txi.reps,tx2gene)

names(txi.genes)

# set the sample names for txi.genes
invisible(sapply(seq(1,3), function(i) {colnames(txi.genes[[i]])<<-mysamples}))

# Read sample metadata
# Data is unordered as it is read in. This means data must be set into the same
# order as the samples were read into mysamples before integrating metadata and
# and read counts

unorderedColData <- read.table("alignment/salmon/Fvenenatum_CarbonRNAseq/F.venenatum/WT_minion/DeSeq2/FvenCarbon_RNAseq_design3.txt",header=T,sep="\t")
colData <- data.frame(unorderedColData[ order(unorderedColData$Sample),])

# Group column
colData$Group <- paste0(colData$Condition,'_', colData$Timepoint)

# Group column
colData$Group2 <- paste0(colData$Condition.1,'_', colData$Timepoint)

# Define the DESeq 'GLM' model
design <- ~ Group
dds <- DESeqDataSetFromTximport(txi.genes,colData,design)

keep <- rowSums(counts(dds)) >= 50
dds <- dds[keep,]

# Library normalisation
dds <- estimateSizeFactors(dds)

# Deseq
dds<-DESeq(dds)

resultsNames(dds)
###
 [1] "Intercept"            "Group_C1_T1_vs_C0_T0" "Group_C1_T2_vs_C0_T0"
 [4] "Group_C1_T3_vs_C0_T0" "Group_C1_T4_vs_C0_T0" "Group_C1_T5_vs_C0_T0"
 [7] "Group_C1_T6_vs_C0_T0" "Group_C1_T7_vs_C0_T0" "Group_C2_T1_vs_C0_T0"
[10] "Group_C2_T2_vs_C0_T0" "Group_C2_T3_vs_C0_T0" "Group_C2_T4_vs_C0_T0"
[13] "Group_C2_T5_vs_C0_T0" "Group_C2_T6_vs_C0_T0" "Group_C2_T7_vs_C0_T0"
[16] "Group_C3_T1_vs_C0_T0" "Group_C3_T2_vs_C0_T0" "Group_C3_T3_vs_C0_T0"
[19] "Group_C3_T4_vs_C0_T0" "Group_C3_T5_vs_C0_T0" "Group_C3_T6_vs_C0_T0"
[22] "Group_C3_T7_vs_C0_T0" "Group_C4_T1_vs_C0_T0" "Group_C4_T2_vs_C0_T0"
[25] "Group_C4_T3_vs_C0_T0" "Group_C4_T4_vs_C0_T0" "Group_C4_T5_vs_C0_T0"
[28] "Group_C4_T6_vs_C0_T0" "Group_C4_T7_vs_C0_T0"
###


# Exploring and exporting results

# res <- results(dds)
# res
# summary(res)

# alpha <- 0.05

# res= results(dds, alpha=alpha,contrast=c("Group","C1_T1","C0_T0"))
# sig.res <- subset(res,padj<=alpha)
# sig.res <- sig.res[order(sig.res$padj),]
# sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
# sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
# summary(sig.res)
# ###
# out of 4192 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 2425, 58%
# LFC < 0 (down)     : 1767, 42%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%
# (mean count < 0)
# ###
# write.table(sig.res,"alignment/salmon/Fvenenatum_CarbonRNAseq_CORRECTED/F.venenatum/WT_minion/corrected/DeSeq2/contrast/C1_vs_C0.txt",sep="\t",na="",quote=F)
# write.table(sig.res.upregulated,"alignment/salmon/Fvenenatum_CarbonRNAseq_CORRECTED/F.venenatum/WT_minion/corrected/DeSeq2/contrast/C1_vs_C0_up.txt",sep="\t",na="",quote=F)
# write.table(sig.res.downregulated,"alignment/salmon/Fvenenatum_CarbonRNAseq_CORRECTED/F.venenatum/WT_minion/corrected/DeSeq2/contrast/C1_vs_C0_down.txt",sep="\t",na="",quote=F)

# res= results(dds, alpha=alpha,contrast=c("Group","C1_T2","C0_T0"))
# sig.res <- subset(res,padj<=alpha)
# sig.res <- sig.res[order(sig.res$padj),]
# sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
# sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
# summary(sig.res)
# ###
# out of 4192 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 2425, 58%
# LFC < 0 (down)     : 1767, 42%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%
# (mean count < 0)
# ###
# write.table(sig.res,"alignment/salmon/Fvenenatum_CarbonRNAseq_CORRECTED/F.venenatum/WT_minion/corrected/DeSeq2/contrast/C2_vs_C0.txt",sep="\t",na="",quote=F)
# write.table(sig.res.upregulated,"alignment/salmon/Fvenenatum_CarbonRNAseq_CORRECTED/F.venenatum/WT_minion/corrected/DeSeq2/contrast/C2_vs_C0_up.txt",sep="\t",na="",quote=F)
# write.table(sig.res.downregulated,"alignment/salmon/Fvenenatum_CarbonRNAseq_CORRECTED/F.venenatum/WT_minion/corrected/DeSeq2/contrast/C2_vs_C0_down.txt",sep="\t",na="",quote=F)

# res= results(dds, alpha=alpha,contrast=c("Condition","C3","C0"))
# sig.res <- subset(res,padj<=alpha)
# sig.res <- sig.res[order(sig.res$padj),]
# sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
# sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
# summary(sig.res)
# ###
# out of 4192 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 2425, 58%
# LFC < 0 (down)     : 1767, 42%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%
# (mean count < 0)
# ###
# write.table(sig.res,"C3_vs_C0.txt",sep="\t",na="",quote=F)
# write.table(sig.res.upregulated,"C3_vs_C0_up.txt",sep="\t",na="",quote=F)
# write.table(sig.res.downregulated,"C3_vs_C0_down.txt",sep="\t",na="",quote=F)

# res= results(dds, alpha=alpha,contrast=c("Condition","C4","C0"))
# sig.res <- subset(res,padj<=alpha)
# sig.res <- sig.res[order(sig.res$padj),]
# sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
# sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
# summary(sig.res)
# ###
# out of 4192 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 2425, 58%
# LFC < 0 (down)     : 1767, 42%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%
# (mean count < 0)
# ###
# write.table(sig.res,"C4_vs_C0.txt",sep="\t",na="",quote=F)
# write.table(sig.res.upregulated,"C4_vs_C0_up.txt",sep="\t",na="",quote=F)
# write.table(sig.res.downregulated,"C4_vs_C0_down.txt",sep="\t",na="",quote=F)


# Sample Distances

# These two are the same
vst<-varianceStabilizingTransformation(dds)
write.csv(assay(vst), file="vst1.csv")
vst5<-varianceStabilizingTransformation(dds,blind=TRUE)
write.csv(assay(vst5), file="vst5.csv")

vst4<-varianceStabilizingTransformation(dds,blind=FALSE)
write.csv(assay(vst4), file="vst4.csv")

vst2<-vst(dds,blind=FALSE)
write.csv(assay(vst2), file="vst2.csv")
vst3<-vst(dds,blind=TRUE)
write.csv(assay(vst3), file="vst3.csv")


pdf("alignment/salmon/Fvenenatum_CarbonRNAseq_CORRECTED/F.venenatum/WT_minion/corrected/heatmap_vst1.pdf", width=12,height=12)
sampleDists<-dist(t(assay(vst)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vst$Condition)
colnames(sampleDistMatrix) <- paste(vst$Timepoint)
colours <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap( sampleDistMatrix, trace="none", col=colours, margins=c(12,12),srtCol=45)
dev.off()

pdf("alignment/salmon/Fvenenatum_CarbonRNAseq_CORRECTED/F.venenatum/WT_minion/corrected/heatmap_vst2.pdf", width=12,height=12)
sampleDists<-dist(t(assay(vst2)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vst$Condition)
colnames(sampleDistMatrix) <- paste(vst$Timepoint)
colours <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap( sampleDistMatrix, trace="none", col=colours, margins=c(12,12),srtCol=45)
dev.off()

pdf("alignment/salmon/Fvenenatum_CarbonRNAseq_CORRECTED/F.venenatum/WT_minion/corrected/heatmap_vst3.pdf", width=12,height=12)
sampleDists<-dist(t(assay(vst3)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vst$Condition)
colnames(sampleDistMatrix) <- paste(vst$Timepoint)
colours <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap( sampleDistMatrix, trace="none", col=colours, margins=c(12,12),srtCol=45)
dev.off()

pdf("alignment/salmon/Fvenenatum_CarbonRNAseq_CORRECTED/F.venenatum/WT_minion/corrected/heatmap_vst4.pdf", width=12,height=12)
sampleDists<-dist(t(assay(vst4)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vst$Condition)
colnames(sampleDistMatrix) <- paste(vst$Timepoint)
colours <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap( sampleDistMatrix, trace="none", col=colours, margins=c(12,12),srtCol=45)
dev.off()

pdf("alignment/salmon/Fvenenatum_CarbonRNAseq_CORRECTED/F.venenatum/WT_minion/corrected/heatmap_vst5.pdf", width=12,height=12)
sampleDists<-dist(t(assay(vst5)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vst$Condition)
colnames(sampleDistMatrix) <- paste(vst$Timepoint)
colours <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap( sampleDistMatrix, trace="none", col=colours, margins=c(12,12),srtCol=45)
dev.off()


# Sample distances measured with rlog transformation:
rld <- rlog(dds)
# pdf("alignment/salmon/Fvenenatum_CarbonRNAseq/F.venenatum/WT_minion/DeSeq2/heatmap_rld.pdf")
# sampleDists <- dist(t(assay(rld)))
# sampleDistMatrix <- as.matrix( sampleDists )
# rownames(sampleDistMatrix) <- paste(rld$Cultivar)
# colnames(sampleDistMatrix) <- paste(rld$Timepoint)
# colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
# heatmap( sampleDistMatrix, trace="none", col=colours, margins=c(12,12),srtCol=45)
# dev.off()

# # MA-plot
# pdf("alignment/salmon/N.ditissima/Hg199/DeSeq2/plotMA_vst.pdf")
# plotMA(res, ylim=c(-2,2))
# dev.off()

# Plot counts
# pdf("alignment/salmon/N.ditissima/Hg199/DeSeq2/plotcounts_dds.pdf")
# plotCounts(dds, gene=which.min(res$padj), intgroup="Cultivar")
# dev.off()

# pdf("alignment/salmon/N.ditissima/Hg199/DeSeq2/plotcounts2_dds.pdf")
# plotCounts(dds, gene=which.min(res$padj), intgroup=c("Cultivar","Timepoint"))
# dev.off()

# PCA plots
pdf("alignment/salmon/Fvenenatum_CarbonRNAseq_CORRECTED/F.venenatum/WT_minion/corrected/PCA_vst5.pdf")
plotPCA(vst5,intgroup=c("Group"))
dev.off()

pdf("alignment/salmon/Fvenenatum_CarbonRNAseq_CORRECTED/F.venenatum/WT_minion/corrected/PCA_vst4.pdf")
plotPCA(vst4,intgroup=c("Group"))
dev.off()

pdf("alignment/salmon/Fvenenatum_CarbonRNAseq_CORRECTED/F.venenatum/WT_minion/corrected/PCA_vst3.pdf")
plotPCA(vst3,intgroup=c("Group"))
dev.off()

pdf("alignment/salmon/Fvenenatum_CarbonRNAseq_CORRECTED/F.venenatum/WT_minion/corrected/PCA_vst2.pdf")
plotPCA(vst2,intgroup=c("Group"))
dev.off()

pdf("alignment/salmon/Fvenenatum_CarbonRNAseq_CORRECTED/F.venenatum/WT_minion/corrected/PCA_vst1.pdf")
plotPCA(vst,intgroup=c("Group"))
dev.off()

# #Plot using rlog transformation:
# pdf("alignment/salmon/Fvenenatum_CarbonRNAseq/F.venenatum/WT_minion/DeSeq2/PCA_rld_group.pdf")
# plotPCA(rld,intgroup=c("Timepoint"))
# dev.off()

#Plot using rlog transformation, showing sample names:

data <- plotPCA(vst, intgroup=c("Condition.1","Timepoint"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
pca_plot<- ggplot(data, aes(PC1, PC2, color=Timepoint, shape=Condition.1)) +
geom_point(size=3) +
xlab(paste0("PC1: ",percentVar[1],"% variance")) +
ylab(paste0("PC2: ",percentVar[2],"% variance")) +
geom_text_repel(aes(label=colnames(vst))) + theme(panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), panel.background = element_blank(),
panel.border = element_rect(colour = "black", fill = NA, size = 1),
axis.text = element_text(size = 14), axis.title = element_text(size = 18))
coord_fixed()
ggsave("PCA_sample_names1.pdf", pca_plot, dpi=300, height=10, width=12)

data <- plotPCA(vst4, intgroup=c("Condition.1","Timepoint"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
pca_plot<- ggplot(data, aes(PC1, PC2, color=Timepoint, shape=Condition.1)) +
geom_point(size=3) +
xlab(paste0("PC1: ",percentVar[1],"% variance")) +
ylab(paste0("PC2: ",percentVar[2],"% variance")) +
geom_text_repel(aes(label=colnames(vst4))) + theme(panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), panel.background = element_blank(),
panel.border = element_rect(colour = "black", fill = NA, size = 1),
axis.text = element_text(size = 14), axis.title = element_text(size = 18))
coord_fixed()
ggsave("PCA_sample_names4.pdf", pca_plot, dpi=300, height=10, width=12)



# data <- plotPCA(vst, intgroup=c("Condition.1","Timepoint"), returnData=TRUE)
# percentVar <- round(100 * attr(data, "percentVar"))
# pca_plot<- ggplot(data, aes(PC1, PC2, color=Timepoint)) +
# geom_point(size=3) +
# xlab(paste0("PC1: ",percentVar[1],"% variance")) +
# ylab(paste0("PC2: ",percentVar[2],"% variance")) + geom_text_repel(aes(label=colnames(vst)))
# coord_fixed()
# ggsave("PCA_sample_names5.pdf", pca_plot, dpi=300, height=10, width=12)


# data <- plotPCA(rld, intgroup=c("Condition"), returnData=TRUE)
# percentVar <- round(100 * attr(data, "percentVar"))
# pca_plot<- ggplot(data, aes(PC1, PC2, color=Condition)) +
# geom_point(size=3) +
# xlab(paste0("PC1: ",percentVar[1],"% variance")) +
# ylab(paste0("PC2: ",percentVar[2],"% variance")) + geom_text_repel(aes(label=colnames(rld)))
# coord_fixed()
# ggsave("PCA_rld_sample_names.pdf", pca_plot, dpi=300, height=10, width=12)


# Gene clustering plots

topVarGenes <-head(order(rowVars(assay(rld)),decreasing=TRUE),50)
pdf("alignment/salmon/N.ditissima/Hg199/DeSeq2/heatmap_new.pdf")
heatmap.2(assay(rld)[topVarGenes,],ColSideColors=c("grey","dodgerblue")[ rld$Timepoint ],scale='row',trace="none",dendrogram="column",col=colorRampPalette(rev(brewer.pal(9,"RdBu")))(255))
dev.off()

topVarGenes <-head(order(rowVars(assay(rld)),decreasing=TRUE),50)
pdf("alignment/salmon/N.ditissima/Hg199/DeSeq2/heatmap_new2.pdf")
colors <- colorRampPalette( rev(brewer.pal(9, "PuOr")) )(255)
sidecols <- c("grey","dodgerblue")[ rld$Cultivar ]
mat <- assay(rld)[ topVarGenes, ]
#mat <- mat - rowMeans(mat)
#colnames(mat) <- paste0(rld$Cultivar,"-",rld$Timepoint)
heatmap.2(mat, trace="none", col=colors, dendrogram="column",ColSideColors=sidecols,labRow=TRUE, mar=c(10,2), scale="row")
dev.off()
```

```r
# Done but not needed
# Make a table of raw counts, normalised counts and fpkm values:
raw_counts <- data.frame(counts(dds, normalized=FALSE))
colnames(raw_counts) <- paste(colData$Sample)
write.table(raw_counts,"alignment/salmon/Fvenenatum_CarbonRNAseq/F.venenatum/WT_minion/DeSeq2/raw_counts.txt",sep="\t",na="",quote=F)
norm_counts <- data.frame(counts(dds, normalized=TRUE))
colnames(norm_counts) <- paste(colData$Sample)
write.table(norm_counts,"alignment/salmon/Fvenenatum_CarbonRNAseq/F.venenatum/WT_minion/DeSeq2/normalised_counts.txt",sep="\t",na="",quote=F)

# robust may be better set at false to normalise based on total counts rather than 'library normalisation factors'
fpkm_counts <- data.frame(fpkm(dds, robust = TRUE))
colnames(fpkm_counts) <- paste(colData$Sample)
write.table(fpkm_counts,"alignment/salmon/Fvenenatum_CarbonRNAseq/F.venenatum/WT_minion/DeSeq2/fpkm_norm_counts.txt",sep="\t",na="",quote=F)
fpkm_counts <- data.frame(fpkm(dds, robust = FALSE))
colnames(fpkm_counts) <- paste(colData$Sample)
write.table(fpkm_counts,"alignment/salmon/Fvenenatum_CarbonRNAseq/F.venenatum/WT_minion/DeSeq2/fpkm_counts.txt",sep="\t",na="",quote=F)


pca(experiment.table="raw_counts.txt", type="counts",
      legend.position="topleft", covariatesInNames=FALSE, samplesName=TRUE,
      principal.components=c(1,2), pdf = TRUE,
      output.folder=getwd())


pca(experiment.table="normalised_counts.txt", type="counts",
      legend.position="topleft", covariatesInNames=FALSE, samplesName=TRUE,
      principal.components=c(1,2), pdf = TRUE,
      output.folder=getwd())

pca(experiment.table="fpkm_counts.txt", type="FPKM",
      legend.position="topleft", covariatesInNames=FALSE, samplesName=TRUE,
      principal.components=c(1,2), pdf = TRUE,
      output.folder=getwd())

pca(experiment.table="fpkm_norm_counts.txt", type="FPKM",
      legend.position="topleft", covariatesInNames=FALSE, samplesName=TRUE,
      principal.components=c(1,2), pdf = TRUE,
      output.folder=getwd())


write.csv(vst, file="vst_all.csv")
write.csv(assay(vst), file="vst_all.csv")
```

```r
# For cosistently expressed genes

topVarGenes <- head( order( rowVars( assay(vst) ), decreasing=FALSE ), 1000 ) # decreasing=TRUE for most variable genes
mat <- assay(vst)[ topVarGenes, ]
write.table(mat,"Topvar2.txt",sep="\t",na="",quote=F)

# Combine both datasets and look shared genes
T1<-read.table("Topvar1.txt",header=T,sep="\t")
T2<-read.table("Topvar2.txt",header=T,sep="\t")
T3<-merge(T1,T2, by.x="ID",by.y="ID",all.x=TRUE,all.y=TRUE) # Print all
T3<-merge(T1,T2, by.x="ID",by.y="ID",all.x=FALSE,all.y=FALSE) # Shared only
write.table(T3,"common.txt",sep="\t",na="",quote=F)
```

### Analysis of conditions

```r
# Load data from SALMON quasi mapping

# import transcript to gene mapping info
tx2gene <- read.table("alignment/salmon/Fvenenatum_CarbonRNAseq_CORRECTED/F.venenatum/WT_minion/corrected/DeSeq2/trans2gene.txt",header=T,sep="\t")

# import quantification files

txi.reps <- tximport(paste(list.dirs("alignment/salmon/Fvenenatum_CarbonRNAseq_CORRECTED/F.venenatum/WT_minion/corrected/DeSeq2", full.names=T,recursive=F),"/quant.sf",sep=""),type="salmon",tx2gene=tx2gene,txOut=T)

# get the sample names from the folders

mysamples <- list.dirs("alignment/salmon/Fvenenatum_CarbonRNAseq_CORRECTED/F.venenatum/WT_minion/corrected/DeSeq2",full.names=F,recursive=F)

# summarise to gene level. This can be done in the tximport step (txOut=F), but is easier to understand in two steps.
txi.genes <- summarizeToGene(txi.reps,tx2gene)

names(txi.genes)

# set the sample names for txi.genes
invisible(sapply(seq(1,3), function(i) {colnames(txi.genes[[i]])<<-mysamples}))

# Read sample metadata
# Data is unordered as it is read in. This means data must be set into the same
# order as the samples were read into mysamples before integrating metadata and
# and read counts

unorderedColData <- read.table("alignment/salmon/Fvenenatum_CarbonRNAseq/F.venenatum/WT_minion/DeSeq2/FvenCarbon_RNAseq_design3.txt",header=T,sep="\t")
colData <- data.frame(unorderedColData[ order(unorderedColData$Sample),])

# Group column
colData$Group <- paste0(colData$Condition,'_', colData$Timepoint)

# Group column
colData$Group2 <- paste0(colData$Condition.1,'_', colData$Timepoint)

# Define the DESeq 'GLM' model
design <- ~ Condition.1
dds <- DESeqDataSetFromTximport(txi.reps,colData,design)

keep <- rowSums(counts(dds)) >= 50
dds <- dds[keep,]

# Library normalisation
dds <- estimateSizeFactors(dds)

# Deseq
dds<-DESeq(dds)

resultsNames(dds)

[1] "Intercept"                           "Condition.1_Glucose_High_vs_Control"
[3] "Condition.1_Glucose_Low_vs_Control"  "Condition.1_Sucrose_High_vs_Control"
[5] "Condition.1_Sucrose_Low_vs_Control" 

# Exploring and exporting results

res <- results(dds)
res
summary(res)

alpha <- 0.05

res= results(dds, alpha=alpha,contrast=c("Condition.1","Glucose_Low","Control"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
summary(sig.res)
###
# out of 2960 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 2026, 68%
# LFC < 0 (down)     : 934, 32%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%
# (mean count < 1)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results
# ###
write.table(sig.res,"alignment/salmon/Fvenenatum_CarbonRNAseq_CORRECTED/F.venenatum/WT_minion/corrected/DeSeq2/contrast/Glucose_Low_vs_Control.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"alignment/salmon/Fvenenatum_CarbonRNAseq_CORRECTED/F.venenatum/WT_minion/corrected/DeSeq2/contrast/Glucose_Low_vs_Control_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"alignment/salmon/Fvenenatum_CarbonRNAseq_CORRECTED/F.venenatum/WT_minion/corrected/DeSeq2/contrast/Glucose_Low_vs_Control_down.txt",sep="\t",na="",quote=F)

res= results(dds, alpha=alpha,contrast=c("Condition.1","Glucose_High","Control"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
summary(sig.res)
###
# out of 2699 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 1751, 65%
# LFC < 0 (down)     : 948, 35%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%
# (mean count < 1)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results
# ###
write.table(sig.res,"alignment/salmon/Fvenenatum_CarbonRNAseq_CORRECTED/F.venenatum/WT_minion/corrected/DeSeq2/contrast/Glucose_High_vs_Control.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"alignment/salmon/Fvenenatum_CarbonRNAseq_CORRECTED/F.venenatum/WT_minion/corrected/DeSeq2/contrast/Glucose_High_vs_Control_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"alignment/salmon/Fvenenatum_CarbonRNAseq_CORRECTED/F.venenatum/WT_minion/corrected/DeSeq2/contrast/Glucose_High_vs_Control_down.txt",sep="\t",na="",quote=F)

res= results(dds, alpha=alpha,contrast=c("Condition.1","Sucrose_High","Control"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
summary(sig.res)
###
# out of 2576 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 1675, 65%
# LFC < 0 (down)     : 901, 35%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%
# (mean count < 1)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results
# ###
write.table(sig.res,"alignment/salmon/Fvenenatum_CarbonRNAseq_CORRECTED/F.venenatum/WT_minion/corrected/DeSeq2/contrast/Sucrose_High_vs_Control.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"alignment/salmon/Fvenenatum_CarbonRNAseq_CORRECTED/F.venenatum/WT_minion/corrected/DeSeq2/contrast/Sucrose_High_vs_Control_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"alignment/salmon/Fvenenatum_CarbonRNAseq_CORRECTED/F.venenatum/WT_minion/corrected/DeSeq2/contrast/Sucrose_High_vs_Control_down.txt",sep="\t",na="",quote=F)


res= results(dds, alpha=alpha,contrast=c("Condition.1","Sucrose_Low","Control"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
summary(sig.res)
###
# out of 2843 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 1852, 65%
# LFC < 0 (down)     : 991, 35%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%
# (mean count < 1)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results
# ###
write.table(sig.res,"alignment/salmon/Fvenenatum_CarbonRNAseq_CORRECTED/F.venenatum/WT_minion/corrected/DeSeq2/contrast/Sucrose_Low_vs_Control.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"alignment/salmon/Fvenenatum_CarbonRNAseq_CORRECTED/F.venenatum/WT_minion/corrected/DeSeq2/contrast/Sucrose_Low_vs_Control_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"alignment/salmon/Fvenenatum_CarbonRNAseq_CORRECTED/F.venenatum/WT_minion/corrected/DeSeq2/contrast/Sucrose_Low_vs_Control_down.txt",sep="\t",na="",quote=F)

### Glucose_High as control

res= results(dds, alpha=alpha,contrast=c("Condition.1","Glucose_Low","Glucose_High"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
summary(sig.res)
###
# out of 4203 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 2333, 56%
# LFC < 0 (down)     : 1870, 44%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%
# (mean count < 0)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results
# ###
write.table(sig.res,"alignment/salmon/Fvenenatum_CarbonRNAseq_CORRECTED/F.venenatum/WT_minion/corrected/DeSeq2/contrast/Glucose_Low_vs_Glucose_High.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"alignment/salmon/Fvenenatum_CarbonRNAseq_CORRECTED/F.venenatum/WT_minion/corrected/DeSeq2/contrast/Glucose_Low_vs_Glucose_High_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"alignment/salmon/Fvenenatum_CarbonRNAseq_CORRECTED/F.venenatum/WT_minion/corrected/DeSeq2/contrast/Glucose_Low_vs_Glucose_High_down.txt",sep="\t",na="",quote=F)


res= results(dds, alpha=alpha,contrast=c("Condition.1","Sucrose_High","Glucose_High"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
summary(sig.res)
###
# out of 169 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 140, 83%
# LFC < 0 (down)     : 29, 17%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%
# (mean count < 0)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results
# ###
write.table(sig.res,"alignment/salmon/Fvenenatum_CarbonRNAseq_CORRECTED/F.venenatum/WT_minion/corrected/DeSeq2/contrast/Sucrose_High_vs_Glucose_High.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"alignment/salmon/Fvenenatum_CarbonRNAseq_CORRECTED/F.venenatum/WT_minion/corrected/DeSeq2/contrast/Sucrose_High_vs_Glucose_High_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"alignment/salmon/Fvenenatum_CarbonRNAseq_CORRECTED/F.venenatum/WT_minion/corrected/DeSeq2/contrast/Sucrose_High_vs_Glucose_High_down.txt",sep="\t",na="",quote=F)


res= results(dds, alpha=alpha,contrast=c("Condition.1","Sucrose_Low","Glucose_High"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
summary(sig.res)
###
# out of 4609 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 2272, 49%
# LFC < 0 (down)     : 2337, 51%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%
# (mean count < 0)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results
# ###
write.table(sig.res,"alignment/salmon/Fvenenatum_CarbonRNAseq_CORRECTED/F.venenatum/WT_minion/corrected/DeSeq2/contrast/Sucrose_Low_vs_Glucose_High.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"alignment/salmon/Fvenenatum_CarbonRNAseq_CORRECTED/F.venenatum/WT_minion/corrected/DeSeq2/contrast/Sucrose_Low_vs_Glucose_High_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"alignment/salmon/Fvenenatum_CarbonRNAseq_CORRECTED/F.venenatum/WT_minion/corrected/DeSeq2/contrast/Sucrose_Low_vs_Glucose_High_down.txt",sep="\t",na="",quote=F)

res= results(dds, alpha=alpha,contrast=c("Condition.1","Sucrose_Low","Sucrose_High"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
sig.res.upregulated <- sig.res[sig.res$log2FoldChange >=1, ]
sig.res.downregulated <- sig.res[sig.res$log2FoldChange <=-1, ]
summary(sig.res)
###
# out of 4192 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 2150, 51%
# LFC < 0 (down)     : 2042, 49%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%
# (mean count < 0)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results
# ###
write.table(sig.res,"alignment/salmon/Fvenenatum_CarbonRNAseq_CORRECTED/F.venenatum/WT_minion/corrected/DeSeq2/contrast/Sucrose_Low_vs_Sucrose_High.txt",sep="\t",na="",quote=F)
write.table(sig.res.upregulated,"alignment/salmon/Fvenenatum_CarbonRNAseq_CORRECTED/F.venenatum/WT_minion/corrected/DeSeq2/contrast//Sucrose_Low_vs_Sucrose_High_up.txt",sep="\t",na="",quote=F)
write.table(sig.res.downregulated,"alignment/salmon/Fvenenatum_CarbonRNAseq_CORRECTED/F.venenatum/WT_minion/corrected/DeSeq2/contrast//Sucrose_Low_vs_Sucrose_High_down.txt",sep="\t",na="",quote=F)
```

# Generating an TSV file with sequencing information

```bash
#Antismash output correction
cat analysis/secondary_metabolites/antismash/F.venenatum/WT_minion_VP/WT_antismash_results_secmet_genes.tsv | sed 's/;//p' | sed 's/;.*//p' | sed 's/Kin.*//p' > analysis/secondary_metabolites/antismash/F.venenatum/WT_minion_VP/WT_antismash_results_secmet_genes_corrected.tsv

for GeneGff in $(ls gene_pred/codingquarry/F.venenatum/WT_minion/final/final_genes_appended_renamed.gff3); do
Strain=WT_minion
Organism=F.venenatum
Assembly=$(ls repeat_masked/F.venenatum/WT_minion/SMARTdenovo/medaka/*_contigs_softmasked_repeatmasker_TPSI_appended.fa)
TFs=$(ls analysis/transcription_factors/F.venenatum/WT_minion/WT_minion_TF_domains.tsv)
InterPro=$(ls gene_pred/interproscan/F.venenatum/WT_minion/WT_minion_interproscan.tsv)
Antismash=$(ls analysis/secondary_metabolites/antismash/F.venenatum/WT_minion_VP/WT_antismash_results_secmet_genes_corrected.tsv)
#Smurf=$(ls analysis/secondary_metabolites/smurf/F.venenatum/WT_minion/WT_minion_smurf_secmet_genes.tsv) # I added cassis genes manually
SwissProt=$(ls gene_pred/swissprot/F.venenatum/WT_minion/swissprot_vJun2020_tophit_parsed.tbl)
Dir1=$(ls -d alignment/salmon/Fvenenatum_CarbonRNAseq_CORRECTED/F.venenatum/WT_minion/corrected/DeSeq2/contrast)
DEG_Files=$(ls \
$Dir1/Glucose_High_vs_Control.txt \
$Dir1/Glucose_Low_vs_Control.txt \
$Dir1/Sucrose_High_vs_Control.txt \
$Dir1/Sucrose_Low_vs_Control.txt \
$Dir1/Glucose_Low_vs_Glucose_High.txt \
$Dir1/Sucrose_High_vs_Glucose_High.txt \
$Dir1/Sucrose_Low_vs_Glucose_High.txt \
$Dir1/Sucrose_Low_vs_Sucrose_High.txt \
| sed -e "s/$/ /g" | tr -d "\n")
OutDir=analysis/annotation_tables_VP/$Organism/$Strain
mkdir -p $OutDir
GeneFasta=$(ls gene_pred/codingquarry/F.venenatum/WT_minion/final/final_genes_appended_renamed.pep.fasta)
ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Annotation_tables
$ProgDir/build_annot_RNAseq.py --gff_format gff3 --gene_gff $GeneGff --gene_fasta $GeneFasta --TFs $TFs --InterPro $InterPro --DEG_files $DEG_Files --Antismash $Antismash --Swissprot $SwissProt > $OutDir/"$Strain"_withDEGs_gene_table.tsv
done
```
