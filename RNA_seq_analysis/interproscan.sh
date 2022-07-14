#script to run interproscan


infile=final_genes_appended_renamed.pep.fasta
outdir=/home/jconnell/interproscan/latest
mkdir -p $outdir
prodir=/home/jconnell/git_repos/emr_repos/Fv_C-variants/RNA_seq_analysis
sbatch $prodir/run_interproscan.sh $infile $outdir



# Create a file with gene and GO term for topGO analysis


OutDir=./
InterProTSV=/home/jconnell/interproscan/final_genes_appended_renamed.pep.fasta.tsv
ProgDir=/home/jconnell
$ProgDir/GO_table.py --interpro $InterProTSV > $OutDir/experiment_all_gene_GO_annots.tsv
cat $OutDir/experiment_all_gene_GO_annots.tsv | sed 's/.t.*//g' > $OutDir/temp1.tsv
cat $OutDir/experiment_all_gene_GO_annots.tsv | cut -f2 > $OutDir/temp2.tsv
paste $OutDir/temp1.tsv $OutDir/temp2.tsv > $OutDir/experiment_all_gene_GO_annots_geneid.tsv
cp $OutDir/experiment_all_gene_GO_annots_geneid.tsv /home/jconnell/out
rm $OutDir/temp1.tsv
rm $OutDir/temp2.tsv



#GO analysis for downreg genes tp1

#BiocManager::install("topGO")
#BiocManager::install("Rgraphviz")
#BiocManager::install("stringr")

library(ggplot2)
library(dplyr)
library(topGO)
library(stringr)
library(Rgraphviz)


setwd("C:/Users/john.connell/Documents/PhD/YearThree/Thesis/RNA_seq_results/gene_ontolgy")

######Prepare file with gene to GO annotation for topGO

#goterms<-read.csv("fv_go_table.csv", header=TRUE)

#genetogo<-goterms[,c(1,18)]
#write.table(genetogo, file = "fv_genetogo.txt", sep = "\t",
#            row.names = FALSE,col.names = FALSE, quote=FALSE)

#GO annotation
gene2GO <- readMappings(file ="fv_gene_go.txt")



#All genes expressed in experiment (genes not expressed in any sample have been removed)
all <-as.data.frame(read.csv("all_tp2.txt",header=TRUE,na.strings="NA"))

# DE genes
de <-as.data.frame(read.csv("all_down_tp2.txt",header=TRUE,na.strings="NA"))

###########################Create TopGO object
interesting <-as.character(de$gene)
allgenes <-as.factor(all$gene)
names(allgenes) <- allgenes
genelist <- factor(as.integer(allgenes %in% interesting))
names(genelist) <- allgenes
GOdata <- new("topGOdata",
              description = "TP2_DR", ontology = "BP",
              allGenes = genelist, annot =annFUN.gene2GO, gene2GO = gene2GO)




###################Run enrichment test
resultFisher <- runTest(GOdata, algorithm="classic", statistic="fisher")
####table of results
sig.tab <- GenTable(GOdata, Fis = resultFisher, topNodes = 15)
sig.tab$Fis<-as.numeric(sig.tab$Fis)
sig.tab$Significant<-as.numeric(sig.tab$Significant)
##Write data of significant GO terms
write.csv(sig.tab, "GOenrichment_GD.csv")
#####Visualising significant terms
sig.tab %>%
  mutate(hitsPerc=(Significant/Annotated)*100) %>%
  ggplot(aes(x=hitsPerc,
             y=Term,
             colour=Fis,
             size=Significant)) +
  geom_point() + theme_bw()+ggtitle("Interaction")+
  theme(axis.text.x=element_text(color="black", size = 15), axis.text.y = element_text(color = "black", size = 15), axis.title.x = element_text(color = "black", size = 25), 
        axis.title.y = element_text(color = "black", size = 25), title = element_text(size = 25, face = "bold")) +
  expand_limits(x=0) +  guides(colour = guide_colourbar(order = 1))+
  labs(x="Hits (%)", y="GO term", colour="p-value", size="Number of genes")

ggsave("GO_analysis_tp2_downreg.jpeg", width = 40, height = 25, units = "cm")
