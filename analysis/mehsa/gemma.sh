
#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=02:00:00
#PBS -V

########Subset Geno and sample files to the list of specified individuals

# on compute node, change directory to 'submission directory':
cd $analysis_dir

analysis_dir=/panfs/panasas01/sscm/qh18484/analysis/mehsa
merged_geno=$analysis_dir/All_chr_genotype.bimbam
temp_pheno=$analysis_dir/temp_pheno.txt

gemma -g $merged_geno -p $temp_pheno -gk 1 -o relatedness_all_geno_cases