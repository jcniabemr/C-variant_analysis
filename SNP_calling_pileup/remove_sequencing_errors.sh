#Remove errors calls from VCF file based on SNP calls in the WT genome 
#To generate error SNPs call SNPs in the WT genome using WT reads

for strain in C2 C3 C4 C5 C6 C7 C8 C9 C10 C11 C12 C13 C14 C15 C16 C17 C18 C19; do
  VCF_file=/projects/fusarium_venenatum_miseq/SNP_calling/pileup_calls/correct_secondary/${strain}_SNP_calls.g.vcf
  Error_call_file=/projects/fusarium_venenatum_miseq/SNP_calling/pileup_calls/correct_secondary/C1_SNP_calls.g.vcf
  OutDir=/projects/fusarium_venenatum_miseq/SNP_calling/pileup_calls/correct_secondary/corrected/${strain}_no_error.vcf  
  ProgDir=/home/connellj/git_repos/emr_repos/Fv_C-variants/SNP_calling_pileup
  $ProgDir/python_remove_sequencing_errors.py --VCF_file $VCF_file --error_SNPs $Error_call_file --outfile $OutDir 
done  



for strain in C1; do
  VCF_file=/projects/fusarium_venenatum_miseq/SNP_calling/pileup_calls/correct_secondary/${strain}_SNP_calls.g.vcf
  Error_call_file=/projects/fusarium_venenatum_miseq/SNP_calling/pileup_calls/correct_secondary/C2_SNP_calls.g.vcf
  OutDir=/projects/fusarium_venenatum_miseq/SNP_calling/pileup_calls/correct_secondary/corrected/${strain}_no_error.vcf  
  ProgDir=/home/connellj/git_repos/emr_repos/Fv_C-variants/SNP_calling_pileup
  $ProgDir/python_remove_sequencing_errors.py --VCF_file $VCF_file --error_SNPs $Error_call_file --outfile $OutDir 
done  



