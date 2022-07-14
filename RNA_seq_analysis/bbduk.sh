#!/usr/bin/env bash
#SBATCH -J bbduk
#SBATCH --partition=long
#SBATCH --mem-per-cpu=12G
#SBATCH --cpus-per-task=12


# Script for bbduk program: Decontamination of rRNA reads in RNAseq data


RIBOKMERS=$1
FORWARD=$2 
REVERSE=$3 
STRAIN=$4
OUTDIR=$5


WorkDir=/home/jconnell/projects/niab/johnc/${SLURM_JOB_USER}_${SLURM_JOBID}
mkdir -p $WorkDir


cp $RIBOKMERS $WorkDir
cp $FORWARD $WorkDir
cp $REVERSE $WorkDir
cd $WorkDir



echo "STRAIN" $STRAIN


RRNA=( k=31 t=10 )


echo "FORWARD" $FORWARD
echo "REVERSE" $REVERSE

F=$(sed 's/.*\///' <<<$FORWARD)
R=$(sed 's/.*\///' <<<$REVERSE)

echo "F" $F
echo "R" $R

# Input files

F_IN=$FORWARD
R_IN=$REVERSE

echo "F_IN" $F_IN
echo "R_IN" $R_IN

# Output files

FOUT=$(sed 's/trim.fq.gz//' <<<$F)
ROUT=$(sed 's/trim.fq.gz//' <<<$R)

OUT=F/"$FOUT"cleaned.fq.gz
OUT2=R/"$ROUT"cleaned.fq.gz
OUTM=F/"$FOUT"rRNA.fq.gz
OUTM2=R/"$ROUT"rRNA.fq.gz
STATS=$STRAIN.stats.txt

echo "OUT" $OUT
echo "OUT2" $OUT2
echo "OUTM" $OUTM
echo "OUTM2" $OUTM2
echo "STATS" $STATS


bbduk=/home/jconnell/miniconda3/bin/bbduk.sh

$bbduk threads=10 in=$F_IN in2=$R_IN out=$OUT out2=$OUT2 outm=$OUTM outm2=$OUTM2 ref=$RIBOKMERS ${RRNA[@]} stats=$STATS

echo "Finish and clean"


mkdir -p $OUTDIR/F
mkdir -p $OUTDIR/R 


cp -r $WorkDir/F/*.fq.gz $OUTDIR/F
cp -r $WorkDir/R/*.fq.gz $OUTDIR/R
rm -r $WorkDir





