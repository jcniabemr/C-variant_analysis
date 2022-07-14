#!/usr/bin/env bash
#SBATCH -J FastQ-mcf
#SBATCH --partition=long
#SBATCH --mem-per-cpu=10G
#SBATCH --cpus-per-task=20



########################################################################
#Input 
# 1st argument: Forward read 
# 2nd argument: Reverse read  
# 3rd argument: Illumina adapters
# 4rd argument: Output directory 
#Output
# Will filter poor quality reads, perform trimming and
# remove illumina adapters.
# rna_qc_fastq-mcf <RNASeq_F.fq> <RNASeq_R.fq> <illumina_adapters.fa> [DNA/RNA]


F_read=$1
R_read=$2
illumina_adapters=$3
F_outdir=$4
R_outdir=$5
Cvariant=$6



WorkDir=/projects/fusarium_venenatum_miseq/${SLURM_JOB_USER}_${SLURM_JOBID}
mkdir -p $WorkDir


cp $F_read $WorkDir
cp $R_read $WorkDir
cp $illumina_adapters $WorkDir
cd $WorkDir



gunzip *R1* 
gunzip *R2*


fastq-mcf \
$illumina_adapters \
*R1* \
*R2* \
-o $WorkDir/"$Cvariant"_forward_trimmed.fastq \
-o $WorkDir/"$Cvariant"_reverse_trimmed.fastq \
-C 1000000 \
-u \
-k 20 \
-t 0.01 \
-q 30 \
-p 5 


gzip $WorkDir/"$Cvariant"_forward_trimmed.fastq
gzip $WorkDir/"$Cvariant"_reverse_trimmed.fastq
cp -r $WorkDir/"$Cvariant"_forward_trimmed.fastq.gz $F_outdir
cp -r $WorkDir/"$Cvariant"_reverse_trimmed.fastq.gz $R_outdir
rm -r $WorkDir
