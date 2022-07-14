#!/usr/bin/env bash
#SBATCH -J ensembl-vep
#SBATCH --partition=long 
#SBATCH --mem-per-cpu=12G
#SBATCH --cpus-per-task=30


##########################################################################
#INPUT:
# 1st argument: snp vcf 
#OUTPUT:
# VCF out


vcf=$1
outdir=$2


WorkDir=/projects/fusarium_venenatum_miseq/${SLURM_JOB_USER}_${SLURM_JOBID}
mkdir -p $WorkDir


cp $vcf $WorkDir
cd $WorkDir


vep=/home/connellj/ensembl-vep/vep
$vep \
-i $vcf \
-o $outdir 


cp $WorkDir/* $outdir
rm -r $WorkDir

