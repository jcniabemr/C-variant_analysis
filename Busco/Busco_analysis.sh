#!/usr/bin/env bash


export MYCONDAPATH=/home/connellj/miniconda2
source ${MYCONDAPATH}/bin/activate py38


for strain in WT; do
 for Assembly in $(ls /projects/fusarium_venenatum_miseq/genomes/$strain/WT_contigs_unmasked.fa); do
   ProgDir=/home/connellj/git_repos/emr_repos/Fv_C-variants/Busco
   BuscoDB=$(ls -d /home/connellj/sordariomycetes_odb10/sordariomycetes_odb10)
   OutDir=/projects/fusarium_venenatum_miseq/genomes/$strain/busco
   mkdir -p $OutDir
   sbatch $ProgDir/busco.sh $Assembly $BuscoDB $OutDir $strain
 done
done  


conda deactivate 










#Run Busco on C-variant genome sequnces.
#BUSCO v5 was run by creating a conda environment with the folowing installed:
#Python v3.8, BioPython, pandas, tBLASTn 2.2+, Augustus 3.2, Prodigal, Metaeuk, HMMER3.1+, SEPP, and R + ggplot2. 
#Finally install busco "conda install -c bioconda -c conda-forge busco=5.1.3"


#1.) Run BUSCO on genomes 















##2.) Create a list of all BUSCO IDs

#OutDir=/projects/fusarium_venenatum_miseq/busco/busco_phylogeny
#mkdir -p $OutDir
#BuscoDb=/home/connellj/sordariomycetes_odb10/sordariomycetes_odb10
#ls -1 $BuscoDb/hmms/*hmm | rev | cut -f1 -d '/' | rev | sed -e 's/.hmm//' > $OutDir/all_buscos_"$BuscoDb".txt





# For BUSCO version 4
#screen -a
#srun --partition long --mem-per-cpu 10G --cpus-per-task 10 --pty bash
# Create a folder for each busco gene
#mkdir temp_busco
#printf "" > /projects/fusarium_venenatum_miseq/busco/busco_phylogeny/single_hits.txt
#  for Busco in $(cat /projects/fusarium_venenatum_miseq/busco/busco_phylogeny/all_buscos_*.txt); do
#  echo $Busco
#  OutDir=/projects/fusarium_venenatum_miseq/busco/busco_phylogeny_2/$Busco
#  mkdir -p $OutDir
  # Move all single copy genes to each folder and rename gene headers
#    for Fasta in $(ls gene_pred/busco/$Organism/$Strain/*/*/*/*/single_copy_busco_sequences/$Busco*.fna); do
#      Strain=$(echo $Fasta | rev | cut -f7 -d '/' | rev)
#      Organism=$(echo $Fasta | rev | cut -f8 -d '/' | rev)
#      FileName=$(basename $Fasta)
#      contig=$(cat $Fasta | grep '>' | sed 's/ <unknown description>//g' | sed 's/>//g')
#      echo ">$Busco:$Strain:$contig" > temp_busco/"$Busco"_"$Strain"_new_names.txt
#      ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools
#      python $ProgDir/replace_fasta_records.py -i $Fasta -r temp_busco/"$Busco"_"$Strain"_new_names.txt -o $OutDir/"$Organism"_"$Strain"_"$Busco".fasta
#      rm temp_busco/"$Busco"_"$Strain"_new_names.txt
    #cat $Fasta | sed "s/:.*.fasta:/:"$Organism"_"$Strain":/g" > $OutDir/"$Organism"_"$Strain"_"$Busco".fasta
#    done
  # Create fasta file containing all busco for alignment
#  cat $OutDir/*_*_"$Busco".fasta > $OutDir/"$Busco"_appended.fasta
#  SingleBuscoNum=$(cat $OutDir/"$Busco"_appended.fasta | grep '>' | wc -l)
#  printf "$Busco\t$SingleBuscoNum\n" >> analysis/popgen/busco_phylogeny/single_hits.txt
#  done
#rm -r temp_busco


