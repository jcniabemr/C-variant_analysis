#!/usr/bin/env bash
#SBATCH -J busco
#SBATCH --partition=himem
#SBATCH --mem-per-cpu=6G
#SBATCH --cpus-per-task=20


##########################################################################
#INPUT:
# 1st argument: Assembly fasta   
# 2nd database option 
#OUTPUT:
# Busco analysis 


Assembly=$1
DatabaseOpt=$2
OutDir=$3

### Output folder
Filename=$(basename "$Assembly")
Prefix="${Filename%.*}"


WorkDir=/projects/fusarium_venenatum_miseq/${SLURM_JOB_USER}_${SLURM_JOBID}
mkdir -p $WorkDir
cp $Assembly $WorkDir
cd $WorkDir


busco \
-o $Prefix \
-i $Filename \
-l $DatabaseOpt \
-m geno \
-c 8 \
--augustus_species fusarium_graminearum


rm $Filename
cp -r * $OutDir
rm -r $WorkDir
