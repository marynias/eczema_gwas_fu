##!/bin/bash
HOME=/panfs/panasas01/sscm/qh18484
analysis=$HOME/analysis/colocalization/twas/gtex
scripts=$HOME/bin/eczema_gwas_fu/colocalization/twas
gwas=$HOME/data/gwas/paternoster2015
my_gwas=results.euro.pval.1k
twas=$HOME/bin/fusion_twas-master
utils=$HOME/bin/eczema_gwas_fu/utils
aspupath=$HOME/bin/aSPUpath2
#Download weights for 6 GTEx datasets.
cd $analysis

cd $twas
mkdir WEIGHTS && cd WEIGHTS

#May want to add Esophagus Mucosa in the future.
wget http://gusevlab.org/projects/fusion/weights/GTEx.Whole_Blood.P01.tar.bz2
wget http://gusevlab.org/projects/fusion/weights/GTEx.Spleen.P01.tar.bz2
wget http://gusevlab.org/projects/fusion/weights/GTEx.Skin_Sun_Exposed_Lower_leg.P01.tar.bz2
wget http://gusevlab.org/projects/fusion/weights/GTEx.Skin_Not_Sun_Exposed_Suprapubic.P01.tar.bz2
wget http://gusevlab.org/projects/fusion/weights/GTEx.Cells_Transformed_fibroblasts.P01.tar.bz2
wget http://gusevlab.org/projects/fusion/weights/GTEx.Cells_EBV-transformed_lymphocytes.P01.tar.bz2

for a in *.bz2
do
tar xjf $a
done

#Match marker IDS to the ones provided in the LD reference file.
#First create a merged reference
cat $twas/LDREF/*bim > $twas/LDREF/all.EUR.bim
python $utils/update_rsid.py --bim $twas/LDREF/all.EUR.bim --tab $gwas/results.euro.pval.tsv \
--head Y --chrom 2 --pos 3 --ref 4 --alt 5 >$twas/LDREF/$my_gwas

cd $analysis
#Create GWAS input.
for chrom in {1..22}
do
echo -e "SNP"'\t'"A1"'\t'"A2"'\t'"Z" > ${my_gwas%.1k}.chr${chrom}
done

cat $twas/LDREF/${my_gwas} | awk -v OFS='\t' '{ print $12, $4, $5, $10 >> "results.euro.pval.chr"$2"" }'

##Run TWAS
function run_fusion_assoc {
chrom=$1
tissue=$2
Rscript $twas/FUSION.assoc_test.R \
--sumstats $analysis/${my_gwas%.1k}.chr${chrom} \
--weights $twas/WEIGHTS/${tissue} \
--weights_dir $twas/WEIGHTS/ \
--ref_ld_chr $twas/LDREF/1000G.EUR. \
--chr $chrom \
--out $analysis/${tissue}_${my_gwas%.1k}.twas.chr${chrom}
}

my_tissue=( Whole_Blood.P01.pos Cells_EBV-transformed_lymphocytes.P01.pos Cells_Transformed_fibroblasts.P01.pos Skin_Not_Sun_Exposed_Suprapubic.P01.pos Skin_Sun_Exposed_Lower_leg.P01.pos Spleen.P01.pos )
for i in "${my_tissue[@]}"
do
for chr in {1..22}
do
qsub -v tissue=$i,chrom=$chr $scripts/sub_run_twas_assoc.sh
done
done

#Joint/conditional tests and plots
my_tissue=( Whole_Blood.P01.pos Cells_EBV-transformed_lymphocytes.P01.pos Cells_Transformed_fibroblasts.P01.pos Skin_Not_Sun_Exposed_Suprapubic.P01.pos Skin_Sun_Exposed_Lower_leg.P01.pos Spleen.P01.pos )
for i in "${my_tissue[@]}"
do
for chr in {1..22}
do
qsub -v tissue=$i,chrom=$chr $scripts/sub_run_twas_cond.sh
done
done

#Run coloc
my_tissue=( Whole_Blood.P01.pos Cells_EBV-transformed_lymphocytes.P01.pos Cells_Transformed_fibroblasts.P01.pos Skin_Not_Sun_Exposed_Suprapubic.P01.pos Skin_Sun_Exposed_Lower_leg.P01.pos Spleen.P01.pos )
for i in "${my_tissue[@]}"
do
for chr in {1..22}
do
qsub -v coloc_threshold=0.05,gwasn=103066,tissue=$i,chrom=$chr $scripts/sub_run_twas_coloc.sh
done
done

##Results harvesting
#Collect significant results into one file (pre-conditioning)
head -1 Skin_Sun_Exposed_Lower_leg.P01.pos_results.euro.pval.twas.chr1.top >gtex_twas.all
for a in *top
do
my_lines=$(wc -l $a | cut -d" " -f1)
if [ $my_lines -gt 1 ]
then
tail -n +2 $a >>gtex_twas.all
fi
done

#Collect significant results into one file (after-conditioning)
head -1 Whole_Blood.P01.pos_results.euro.pval.twas.chr20.analysis.joint_included.dat >gtex_twas.independent
for a in *_included.dat
do
tail -n +2 $a >>gtex_twas.independent
done

#Collect significant coloc results
head -1 Skin_Not_Sun_Exposed_Suprapubic.P01.pos_results.euro.pval.coloc.chr22 >coloc_twas.sig0.45
for a in *coloc*
do
cat $a | awk '$25 != "NA" {print $0}' | awk -v OFS="\t" '$25 > 0.45 {print $0}' | tail -n +2 >>coloc_twas.sig0.45
done

