#!/bin/bash
#PBS -l nodes=2:ppn=1
#PBS -l walltime=12:00:00
#PBS -V

cd $PBS_O_WORKDIR

for a in {1..22}
do
plink --bfile ALL.chr${a}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_EUR_no_mono \
--exclude EUR_all_1k-merge.missnp --make-bed --out ALL.chr${a}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_EUR_no_mono_tmp
done

plink --bmerge ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_EUR_no_mono_tmp \
--bmerge ALL.chr2.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_EUR_no_mono_tmp \
--bmerge ALL.chr3.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_EUR_no_mono_tmp \
--bmerge ALL.chr4.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_EUR_no_mono_tmp \
--bmerge ALL.chr5.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_EUR_no_mono_tmp \
--bmerge ALL.chr6.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_EUR_no_mono_tmp \
--bmerge ALL.chr7.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_EUR_no_mono_tmp \
--bmerge ALL.chr8.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_EUR_no_mono_tmp \
--bmerge ALL.chr9.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_EUR_no_mono_tmp \
--bmerge ALL.chr10.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_EUR_no_mono_tmp \
--bmerge ALL.chr11.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_EUR_no_mono_tmp \
--bmerge ALL.chr12.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_EUR_no_mono_tmp \
--bmerge ALL.chr13.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_EUR_no_mono_tmp \
--bmerge ALL.chr14.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_EUR_no_mono_tmp \
--bmerge ALL.chr15.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_EUR_no_mono_tmp \
--bmerge ALL.chr16.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_EUR_no_mono_tmp \
--bmerge ALL.chr17.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_EUR_no_mono_tmp \
--bmerge ALL.chr18.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_EUR_no_mono_tmp \
--bmerge ALL.chr19.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_EUR_no_mono_tmp \
--bmerge ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_EUR_no_mono_tmp \
--bmerge ALL.chr21.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_EUR_no_mono_tmp \
--bmerge ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_EUR_no_mono_tmp \
--make-bed --out merged

plink --merge-list 1k_eur_merge_files_temp.txt --make-bed --out EUR_all_1k
rm -rf ALL.chr$a.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_EUR_no_mono_tmp*