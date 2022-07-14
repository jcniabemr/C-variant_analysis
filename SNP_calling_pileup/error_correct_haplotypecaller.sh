#!/usr/bin/env bash
#SBATCH -J haplotype_caller
#SBATCH --partition=long 
#SBATCH --mem-per-cpu=12G
#SBATCH --cpus-per-task=30



##########################################################################
#INPUT:
# 1st argument: Refrence fasta 
# 2nd argument: sample name (prefix) to be used to identify it in the future
# 3nd argument: Input BAM 
#OUTPUT:
# VCF out


reference=$1
alignment=$2
outdir=$3



WorkDir=/projects/fusarium_venenatum_miseq/${SLURM_JOB_USER}_${SLURM_JOBID}
mkdir -p $WorkDir


cp $reference $WorkDir
cp $alignment $WorkDir
cd $WorkDir


samtools faidx *.fa
samtools index *rg.bam



picard=/home/connellj/miniconda2/share/picard-2.18.29-0/picard.jar
java -jar $picard CreateSequenceDictionary \
	R=WT_contigs_unmasked.fa \
	O=WT_contigs_unmasked.dict 



gatk=/scratch/software/GenomeAnalysisTK-3.6
java -jar $gatk/GenomeAnalysisTK.jar \
     -T HaplotypeCaller \
     -R WT_contigs_unmasked.fa \
     -ploidy 1 \
     -nct 24 \
     --allow_potentially_misencoded_quality_scores \
     -I WT_contigs_unmasked.fa_aligned_sorted_nomulti_proper_sorted_nodup_rg.bam \
     -o error_snps.vcf


cp $WorkDir/*.vcf $outdir
rm -r $WorkDir