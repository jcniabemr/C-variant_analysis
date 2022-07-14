#!/usr/bin/env bash
#SBATCH -J indel_realignment
#SBATCH --partition=short
#SBATCH --mem-per-cpu=12G
#SBATCH --cpus-per-task=30


##########################################################################
#INPUT:
# 1st argument: Refrence fasta 
# 2nd argument: sample name (prefix) to be used to identify it in the future
# 3rd argument: input BAM file from pre_snp_calling file with your mappings with duplicates marked no multimapping sorted 
# 4th argument: target intervals required for indel rearangemenr 
#OUTPUT:
# realigned bam files prefixed with strain ID

reference=$1
strain=$2
input_bam=$3
target_intervals=$4
outdir=$5


WorkDir=/projects/fusarium_venenatum_miseq/${SLURM_JOB_USER}_${SLURM_JOBID}
mkdir -p $WorkDir


cp $reference $WorkDir
cp $input_bam $WorkDir
cp $target_intervals $WorkDir
cd $WorkDir


samtools faidx WT_contigs_unmasked.fa
samtools index *_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam


picard=/home/connellj/miniconda2/share/picard-2.18.29-0/picard.jar
java -jar $picard CreateSequenceDictionary \
	R=WT_contigs_unmasked.fa \
	O=WT_contigs_unmasked.dict 


gatk=/scratch/software/GenomeAnalysisTK-3.6
java -jar $gatk/GenomeAnalysisTK.jar \
     -T IndelRealigner \
     -R WT_contigs_unmasked.fa \
     -I *_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -targetIntervals realigner.intervals \
     -o "$strain"_realigned.bam 


cp $WorkDir/"$strain"_realigned.bam $outdir
cp $WorkDir/realigner.intervals $outdir
rm -r $WorkDir