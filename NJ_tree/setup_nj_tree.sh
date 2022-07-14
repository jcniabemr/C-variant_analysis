#Create a NJ tree using variant call data 
#For my run I had to install the following:
	#1 "conda install -c conda-forge r-ape"
	#2 "conda install -c bioconda perl-vcftools-vcf"


VCF_file=/projects/fusarium_venenatum_miseq/SNP_calling/pileup_calls/WT_contigs_unmaskedSNPs_filtered.vcf
ploidy=1
progdir=/home/connellj/git_repos/emr_repos/Fv_C-variants/NJ_tree
sbatch $progdir/create_nj_tree.sh $VCF_file $ploidy