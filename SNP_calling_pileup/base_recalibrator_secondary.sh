#!/usr/bin/env bash
#SBATCH -J base_recalibrator_seconsary
#SBATCH --partition=short
#SBATCH --mem-per-cpu=12G
#SBATCH --cpus-per-task=30


##########################################################################
#INPUT:
# 1st argument: Refrence fasta 
# 2nd argument: SNP call vcf for variant masking 
# 3rd argument: indell call vcf for variant masking 
# 4th argument: sv call vcf for masking 
# 5th argument: sample name (prefix) to be used to identify it in the future
# 6th argument: input BAM file from pre_snp_calling file with your mappings with duplicates marked no multimapping sorted  
# 7th argument: input recal table generated from primary base recalibation
#OUTPUT:
# Secondary realibration table used later for plotting improved root mean square error 


reference=$1
SNP=$2
INDEL=$3
SV=$4
strain=$5
input_bam=$6
recal_table=$7
outdir=$8


WorkDir=/projects/fusarium_venenatum_miseq/${SLURM_JOB_USER}_${SLURM_JOBID}
mkdir -p $WorkDir


cp $reference $WorkDir
cp $SNP $WorkDir
cp $INDEL $WorkDir
cp $SV $WorkDir
cp $input_bam $WorkDir
cp $recal_table $WorkDir
cd $WorkDir


samtools faidx WT_contigs_unmasked.fa
samtools index *_realigned.bam


picard=/home/connellj/miniconda2/share/picard-2.18.29-0/picard.jar
java -jar $picard CreateSequenceDictionary \
	R=WT_contigs_unmasked.fa \
	O=WT_contigs_unmasked.dict 


gatk=/scratch/software/GenomeAnalysisTK-3.6
java -jar $gatk/GenomeAnalysisTK.jar \
     -T BaseRecalibrator \
     -R WT_contigs_unmasked.fa \
     -I "$strain"_realigned.bam \
     -knownSites corrected_snp.bam \
     -knownSites Fven_svaba_sv.svaba.unfiltered.indel.vcf \
     -knownSites Fven_svaba_sv.svaba.unfiltered.sv.vcf \
     -BQSR "$strain"_recal.table \
     -o "$strain"_secondary_recal.table 


cp $WorkDir/"$strain"_secondary_recal.table $outdir
rm -r $WorkDir