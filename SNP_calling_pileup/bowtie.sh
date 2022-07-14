#!/usr/bin/env bash
#SBATCH -J bowtie
#SBATCH --partition=long
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=8



Assembly=$(basename $1)
Read_F=$(basename $2)
Read_R=$(basename $3)
OutDir=$4



if [ $5 ]; then

  AddOptions="--rg-id $5"

else

  AddOptions=""

fi


CurDir=$PWD


WorkDir=$TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}

mkdir -p $WorkDir

cd $WorkDir

cp $CurDir/$1 $Assembly

cp $CurDir/$2 $Read_F

cp $CurDir/$3 $Read_R


bowtie2-build $Assembly $Assembly.indexed

bowtie2 -p 8 -X 1200 --no-mixed $AddOptions -x $Assembly.indexed  -1 $Read_F -2 $Read_R  -S "$Assembly"_aligned.sam 2>&1 | tee bowtie_log.txt

samtools view --threads 8 -bS "$Assembly"_aligned.sam -o "$Assembly"_aligned.bam

samtools sort --threads 8 -o "$Assembly"_aligned_sorted.bam "$Assembly"_aligned.bam

samtools index -@ 8 "$Assembly"_aligned_sorted.bam "$Assembly"_aligned_sorted.bam.index

samtools view --threads 8 -f 2 -F 104 "$Assembly"_aligned_sorted.bam | cut -f3 | uniq -c > "$Assembly"_RPK.txt



cp -r $WorkDir/WT_contigs_unmasked.fa_aligned_sorted* $OutDir

rm -r $WorkDir