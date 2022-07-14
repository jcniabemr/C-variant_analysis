#This document details commands used for SNP calling pileup in F. venenatum C-variants vs a reference wild-type assembly.


#1.) Carryout alignment of the WT raw sequencing reads to the assembled WT genome using bowtie. 


Reference=../../projects/fusarium_venenatum_miseq/genomes/WT/WT_contigs_unmasked.fa
for StrainPath in $(ls -d ../../projects/fusarium_venenatum_miseq/genomes/WT); do
  F_Read=$(ls $StrainPath/F/trim/FvenWT_S1_L001_R1_001_trim.fq.gz)
  R_Read=$(ls $StrainPath/R/trim/FvenWT_S1_L001_R2_001_trim.fq.gz)
  OutDir=../../projects/fusarium_venenatum_miseq/SNP_calling/error_correction 
  mkdir -p $OutDir
  ProgDir=/home/connellj/git_repos/emr_repos/Fv_C-variants/SNP_calling_pileup
  sbatch $ProgDir/bowtie.sh $Reference $F_Read $R_Read $OutDir
  cd $OutDir 
  samtools view -h -o WT_contigs_unmasked.fa_aligned_sorted.sam WT_contigs_unmasked.fa_aligned_sorted.bam
done 


#2.) Remove multimapping reads, discordant reads. PCR and optical duplicates, and add read group and sample name to each mapped read. 


input=/projects/fusarium_venenatum_miseq/SNP_calling/error_correction/WT_contigs_unmasked.fa_aligned_sorted.sam
OutDir=/projects/fusarium_venenatum_miseq/SNP_calling/error_correction
ProgDir=/home/connellj/git_repos/emr_repos/Fv_C-variants/SNP_calling_pileup
sbatch $ProgDir/pre_SNP_calling_error.sh $input $OutDir 
 


#2.) Submit GATK SNP calling. 

Reference=/projects/fusarium_venenatum_miseq/genomes/WT/WT_contigs_unmasked.fa  
alignment=/projects/fusarium_venenatum_miseq/SNP_calling/error_correction/WT_contigs_unmasked.fa_aligned_sorted_nomulti_proper_sorted_nodup_rg.bam 
OutDir=/projects/fusarium_venenatum_miseq/SNP_calling/error_correction
ProgDir=/home/connellj/git_repos/emr_repos/Fv_C-variants/SNP_calling_pileup
sbatch $ProgDir/error_correct_haplotypecaller.sh $Reference $alignment $OutDir 


#Remove errors calls from VCF file based on SNP calls in the WT genome 
#To generate error SNPs call SNPs in the WT genome using WT reads


  VCF_file=/projects/fusarium_venenatum_miseq/SNP_calling/pileup_calls/WT_contigs_unmaskedSNPs_filtered_nonsyn.vcf
  Error_call_file=/projects/fusarium_venenatum_miseq/SNP_calling/error_correction/error_snps.vcf
  OutDir=/projects/fusarium_venenatum_miseq/SNP_calling/error_correction/no_error.vcf  
  ProgDir=/home/connellj/git_repos/emr_repos/Fv_C-variants/SNP_calling_pileup
  $ProgDir/python_remove_sequencing_errors.py --VCF_file $VCF_file --error_SNPs $Error_call_file --outfile $OutDir 
