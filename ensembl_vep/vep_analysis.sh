#Analyse snp data with vep


vcf=
outdir=/projects/fusarium_venenatum_miseq/vep
progdir=/home/connellj/git_repos/emr_repos/Fv_C-variants/ensembl_vep
sbatch $progdir $vcf $outdir