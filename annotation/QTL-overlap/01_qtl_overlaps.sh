##!/bin/bash
HOME=/panfs/panasas01/sscm/qh18484
analysis=$HOME/analysis/annotation/QTL-overlaps
data_manipulation=$HOME/analysis/annotation/data_manipulation
scripts=$HOME/bin/eczema_gwas_fu/annotation/QTL-overlap
datasets_old=/panfs/panasas01/dedicated-mrcieu/users/msobczyk/Datasets
datasets=$HOME/working/data/Datasets
onek=/panfs/panasas01/sscm/qh18484/analysis/bayesian_fm/RefPanel/1kGenomes
utils=$HOME/bin/eczema_gwas_fu/utils
cd $datasets
#Idenitfy overlaps of SNPs in the 3Mbp interval around index SNP with any SNPs harboring any significant QTL in the downloaded datasets.

#Grep for first rsids, and then for chromosome, position, and alleles.
#Use three sets: 1Mbp around index SNP, 3Mbp around index SNP and r2 >0.2 in 1K Genomes.

#First convert all the excel files into tsv.

for i in $(find . -type f -name \*.xls)
do 
Rscript --vanilla $scripts/convert_xls_to_txt.R $i
done

for i in $(find . -type f -name \*.xlsx)
do 
Rscript --vanilla $scripts/convert_xls_to_txt.R $i
done

Rscript --vanilla $scripts/convert_xls_to_txt.R $datasets/promoter-enhancer/Rinaldi2016/stem_2032_mmc5.xlsb

#Prepare input from Fagny 2017.
cd $datasets/network/Fagny2017
qsub $scripts/sub_fagny2017_extract.sh
cd $datasets
#Merge together all results from pQTL analysis from Suhre et al. (2017)
touch ./pQTL/Suhre2017/all_suhre2017.txt
for a in ./pQTL/Suhre2017/*gz
do
gunzip $a
done

cat ./pQTL/Suhre2017/*out >>./pQTL/Suhre2017/all_suhre2017.txt
#Annotate input from Pala 2017 with 1k rsids.
cd $datasets/eQTL/Pala2017
python $utils/update_rsid_all.py --bim $onek/diagnostics/ALL.all.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_EUR_no_mono.bim --tab eQTLs.Fdr0.05.WithConditional.Annot.Recoded.tsv --head Y --chrom 1 --pos 2 --ref 4 --alt 5 >eQTLs.Fdr0.05.WithConditional.Annot.Recoded.tsv.1k
python $utils/update_rsid_all.py --bim $onek/diagnostics/ALL.all.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_EUR_no_mono.bim --tab isoQTLs.Fdr0.05.WithConditional.Annot.Recoded.tsv --head Y --chrom 1 --pos 2 --ref 4 --alt 5 >isoQTLs.Fdr0.05.WithConditional.Annot.Recoded.tsv.1k

function grep_snp_header {
my_list_file=$1
my_list_name=$2
my_textfile=$3
head -1 $my_textfile > ${my_textfile}_vs_${my_list_name}
grep -w -F -f $my_list_file $my_textfile >>${my_textfile}_vs_${my_list_name}
}

function grep_snp_noheader {
my_list_file=$1
my_list_name=$2
my_textfile=$3
grep -w -F -f $my_list_file $my_textfile >>${my_textfile}_vs_${my_list_name}
}

function grep_snp_header_n {
my_list_file=$1
my_list_name=$2
my_textfile=$3
my_n=$4
head -${my_n} $my_textfile > ${my_textfile}_vs_${my_list_name}
grep -w -F -f $my_list_file $my_textfile >>${my_textfile}_vs_${my_list_name}
}

function grep_chosen_files {
my_list_file=$1
my_list_name=$2
for my_textfile in ./eQTL/Battle2013/Supplemental_Data_1.1.tsv \
./eQTL/Battle2013/Supplemental_Data_1.2.tsv \
./eQTL/Battle2013/Supplemental_Data_1.3.tsv \
./eQTL/Battle2013/Supplemental_Data_1.4.tsv \
./eQTL/Battle2013/Supplemental_Data_1.5.tsv \
./eQTL/Battle2013/Supplemental_Data_1.6.tsv \
./eQTL/Bonder2017/2015_09_02_cis_meQTLsFDR0.05-CpGLevel.txt \
./eQTL/Bonder2017/2015_09_02_trans_meQTLsFDR0.05-CpGLevel.txt \
./eQTL/Bonder2017/exon_level_eQTLs_independent_effects.txt \
./eQTL/Bonder2017/exon-ratio_level_eQTLs_independent_effects.txt \
./eQTL/Bonder2017/gene_level_eQTLs.txt \
./eQTL/Bonder2017/gene_level_eQTLs_independent_effects_interactions.txt \
./eQTL/Bonder2017/ng.3721-S5.1.tsv \
./eQTL/Bonder2017/ng.3721-S5.2.tsv \
./eQTL/Bonder2017/ng.3721-S6.1.tsv \
./eQTL/Bonder2017/ng.3721-S7.1.tsv \
./eQTL/Bonder2017/ng.3721-S8.1.tsv \
./eQTL/Bonder2017/ng.3721-S9.1.tsv \
./eQTL/Bonder2017/polyA-ratio_level_eQTLs_independent_effects.txt \
./eQTL/Brown2017/ng.3979-S3.txt 
do
grep_snp_header $my_list_file $my_list_name $my_textfile
done
}

function grep_chosen_files2 {
my_list_file=$1
my_list_name=$2
for my_textfile in ./eQTL/Chen2016/mono_gene_nor_combat_peer_10_all_summary.txt \
./eQTL/Chen2016/neut_gene_nor_combat_peer_10_all_summary.txt \
./eQTL/Chen2016/tcel_gene_nor_combat_peer_10_all_summary.txt \
./hQTL/Chen2016/mono_K27AC_log2rpm_peer_10_all_summary.txt \
./hQTL/Chen2016/mono_K4ME1_log2rpm_peer_10_all_summary.txt \
./hQTL/Chen2016/neut_K27AC_log2rpm_peer_10_all_summary.txt \
./hQTL/Chen2016/neut_K4ME1_log2rpm_peer_10_all_summary.txt \
./hQTL/Chen2016/tcel_K27AC_log2rpm_peer_10_all_summary.txt \
./hQTL/Chen2016/tcel_K4ME1_log2rpm_peer_10_all_summary.txt 
do
grep_snp_noheader $my_list_file $my_list_name $my_textfile
#Filter significant results.
cat ${my_textfile}_vs_${my_list_name} | awk -v OFS="\t" '$7 < 0.05 {print $0}' >${my_textfile}_vs_${my_list_name}.sig
done
for my_textfile in ./eQTL/Chen2016/tcel_gene_WASP_ASE_all.txt \
./eQTL/Chen2016/neut_gene_WASP_ASE_all.txt \
./eQTL/Chen2016/mono_gene_WASP_ASE_all.txt \
./hQTL/Chen2016/tcel.H3K4me1_peak_WASP_ASE_all.txt \
./hQTL/Chen2016/tcel.H3K27ac_peak_WASP_ASE_all.txt \
./hQTL/Chen2016/neut.H3K4me1_peak_WASP_ASE_all.txt \
./hQTL/Chen2016/neut.H3K27ac_peak_WASP_ASE_all.txt \
./hQTL/Chen2016/mono.H3K4me1_peak_WASP_ASE_all.txt \
./hQTL/Chen2016/mono.H3K27ac_peak_WASP_ASE_all.txt 
do
grep_snp_noheader $my_list_file $my_list_name $my_textfile
#Filter significant results.
cat ${my_textfile}_vs_${my_list_name} | awk -v OFS="\t" '$6 < 0.05 {print $0}' >${my_textfile}_vs_${my_list_name}.sig
done
}

function grep_chosen_files3 {
my_list_file=$1
my_list_name=$2
for my_textfile in ./eQTL/Dimas2009/borel_TableS1_and_S2.2.tsv \
./eQTL/Dimas2009/GENCORD2_EQTL_FDR10_F_183 \
./eQTL/Dimas2009/GENCORD2_EQTL_FDR10_L_185 \
./eQTL/Dimas2009/GENCORD2_EQTL_FDR10_T_186 \
./eQTL/Dimas2009/GENCORD2_MQTL_FDR10_F_107 \
./eQTL/Dimas2009/GENCORD2_MQTL_FDR10_L_111 \
./eQTL/Dimas2009/GENCORD2_MQTL_FDR10_T_66 \
./eQTL/Ding2010/NN57subj_p1E-5_annotate_allcis_1Mb.tbl \
./eQTL/Ding2010/PN53subj_p1E-5_annotate_allcis_1Mb.tbl \
./eQTL/Ding2010/PP53subj_p1E-5_annotate_allcis_1Mb.tbl \
./eQTL/Fairfax2012/ng.2205-S3.1.tsv \
./eQTL/Fairfax2012/ng.2205-S2.1.tsv \
./eQTL/Fairfax2012/ng.2205-S2.2.tsv \
./eQTL/Fairfax2012/ng.2205-S2.3.tsv \
./eQTL/Fairfax2012/ng.2205-S3.2.tsv \
./eQTL/Fairfax2012/ng.2205-S5.1.tsv \
./eQTL/Fairfax2012/ng.2205-S7.1.tsv \
./eQTL/Fairfax2014/1246949stableS2.1.tsv \
./eQTL/Fairfax2014/1246949stableS2.2.tsv \
./eQTL/Fairfax2014/1246949stableS2.3.tsv \
./eQTL/Fairfax2014/1246949stableS3.1.tsv 
do
grep_snp_header $my_list_file $my_list_name $my_textfile
done
}

function grep_chosen_files4 {
my_list_file=$1
my_list_name=$2
for my_textfile in ./eQTL/Kasela2017/journal.pgen.1006643.s013.1.tsv \
./eQTL/Kasela2017/journal.pgen.1006643.s013.2.tsv \
./eQTL/Kim-Hellmuth2017/41467_2017_366_MOESM3_ESM.1.tsv \
./eQTL/Kim-Hellmuth2017/41467_2017_366_MOESM5_ESM.1.tsv \
./eQTL/Lappalainen2013/nature12531-s2.2.tsv \
./eQTL/Lappalainen2013/YRI89.trratio.cis.FDR5.best.rs137.txt \
./eQTL/Lappalainen2013/YRI89.trratio.cis.FDR5.all.rs137.txt \
./eQTL/Lappalainen2013/YRI89.repeat.cis.FDR5.best.rs137.txt \
./eQTL/Lappalainen2013/YRI89.repeat.cis.FDR5.all.rs137.txt \
./eQTL/Lappalainen2013/YRI89.mi.cis.FDR5.best.rs137.txt \
./eQTL/Lappalainen2013/YRI89.mi.cis.FDR5.all.rs137.txt \
./eQTL/Lappalainen2013/YRI89.gene.cis.FDR5.best.rs137.txt \
./eQTL/Lappalainen2013/YRI89.gene.cis.FDR5.all.rs137.txt \
./eQTL/Lappalainen2013/YRI89.exon.cis.FDR5.best.rs137.txt \
./eQTL/Lappalainen2013/YRI89.exon.cis.FDR5.all.rs137.txt \
./eQTL/Lappalainen2013/EUR373.trratio.cis.FDR5.best.rs137.txt \
./eQTL/Lappalainen2013/EUR373.trratio.cis.FDR5.all.rs137.txt \
./eQTL/Lappalainen2013/EUR373.repeat.cis.FDR5.best.rs137.txt \
./eQTL/Lappalainen2013/EUR373.repeat.cis.FDR5.all.rs137.txt \
./eQTL/Lappalainen2013/EUR373.gene.cis.FDR5.best.rs137.txt \
./eQTL/Lappalainen2013/EUR373.gene.cis.FDR5.all.rs137.txt \
./eQTL/Lappalainen2013/EUR373.exon.cis.FDR5.best.rs137.txt \
./eQTL/Lappalainen2013/EUR373.exon.cis.FDR5.all.rs137.txt \
./eQTL/Lappalainen2013/EUR363.mi.cis.FDR5.best.rs137.txt \
./eQTL/Lappalainen2013/EUR363.mi.cis.FDR5.all.rs137.txt \
./eQTL/Momozawa2018/41467_2018_4365_MOESM3_ESM.1.tsv \
./eQTL/Thalayasingam2018/A+R_Supplementary_Tables_as_Submitted.3.tsv \
./eQTL/Thalayasingam2018/A+R_Supplementary_Tables_as_Submitted.4.tsv \
./eQTL/Thalayasingam2018/A+R_Supplementary_Tables_as_Submitted.5.tsv \
./eQTL/Thalayasingam2018/A+R_Supplementary_Tables_as_Submitted.6.tsv \
./eQTL/Walsh2016/13059_2016_948_MOESM1_ESM.txt \
./eQTL/Walsh2016/13059_2016_948_MOESM3_ESM.txt \
./eQTL/Wijst2018/41588_2018_89_MOESM3_ESM.1.tsv \
./eQTL/Wijst2018/41588_2018_89_MOESM3_ESM.2.tsv \
./eQTL/Wijst2018/41588_2018_89_MOESM6_ESM.1.tsv \
./eQTL/Xia2012/eQTL_Qvalue_cutoff_hapmap3_cis_hg19.txt \
./eQTL/Xia2012/eQTL_Qvalue_cutoff_hapmap3_trans_hg19.txt \
./eQTL/Yao2017/1-s2.0-S0002929717300708-mmc2.1.tsv \
./eQTL/Yao2017/1-s2.0-S0002929717300708-mmc3.1.tsv \
./eQTL/Yao2017/1-s2.0-S0002929717300708-mmc4.1.tsv \
./eQTL/Yao2017/1-s2.0-S0002929717300708-mmc5.1.tsv \
./eQTL/Yao2017/1-s2.0-S0002929717300708-mmc6.1.tsv \
./eQTL/Zhang2017/231423-2.13.tsv \
./eQTL/Zhang2017/231423-2.14.tsv \
./eQTL/Zhang2017/231423-2.15.tsv \
./eQTL/Zhang2017/231423-2.16.tsv \
./eQTL/Zhernakova2017/2015_09_02_cis_meQTLsFDR0.05-CpGLevel.txt \
./eQTL/Zhernakova2017/2015_09_02_trans_meQTLsFDR0.05-CpGLevel.txt \
./eQTL/Zhernakova2017/eQTLsFDR0.05-ProbeLevel.txt \
./eQTL/Zhernakova2017/exon_level_eQTLs_independent_effects.txt \
./eQTL/Zhernakova2017/exon-ratio_level_eQTLs_independent_effects.txt \
./eQTL/Zhernakova2017/gene_level_eQTLs.txt \
./eQTL/Zhernakova2017/gene_level_eQTLs_independent_effects_interactions.txt \
./eQTL/Zhernakova2017/polyA-ratio_level_eQTLs_independent_effects.txt \
./hQTL/Pelikan2018/41467_2018_5328_MOESM4_ESM.1.tsv \
./hQTL/Pelikan2018/41467_2018_5328_MOESM6_ESM.1.tsv \
./hQTL/Pelikan2018/41467_2018_5328_MOESM7_ESM.1.tsv \
./hQTL/Pelikan2018/41467_2018_5328_MOESM8_ESM.1.tsv \
./hQTL/Pelikan2018/41467_2018_5328_MOESM9_ESM.1.tsv \
./hQTL/Pelikan2018/41467_2018_5328_MOESM10_ESM.1.tsv \
./promoter-enhancer/Javierre2016/1-s2.0-S0092867416313228-mmc2.2.tsv \
./promoter-enhancer/Javierre2016/1-s2.0-S0092867416313228-mmc2.3.tsv \
./promoter-enhancer/Javierre2016/1-s2.0-S0092867416313228-mmc2.4.tsv \
./promoter-enhancer/Javierre2016/1-s2.0-S0092867416313228-mmc3.7.tsv 
do
grep_snp_header $my_list_file $my_list_name $my_textfile
done
}

function grep_chosen_files5 {
my_list_file=$1
my_list_name=$2
for my_textfile in ./eQTL/Kasela2017/journal.pgen.1006643.s011.1.tsv \
./eQTL/Kasela2017/journal.pgen.1006643.s010.2.tsv \
./eQTL/Kasela2017/journal.pgen.1006643.s010.1.tsv \
./eQTL/Wijst2018/41588_2018_89_MOESM4_ESM.1.tsv \
./eQTL/Wijst2018/41588_2018_89_MOESM5_ESM.1.tsv \
./eQTL/Wijst2018/41588_2018_89_MOESM6_ESM.2.tsv 
do
grep_snp_header_n $my_list_file $my_list_name $my_textfile 2
done
}

function grep_chosen_files6 {
my_list_file=$1
my_list_name=$2
for my_textfile in ./pQTL/Suhre2017/all_suhre2017.txt
do
grep_snp_noheader $my_list_file $my_list_name $my_textfile 
done
}

function grep_godmc {
my_list_file=$1
my_list_name=$2
tail -n +2 ./mQTL/GoDMC/snps_36cohorts16_rsid.txt | sort >./mQTL/GoDMC/snps_36cohorts16_rsid_sorted.txt 
tail -n +2 $my_list_file  | sort >./mQTL/GoDMC/${my_list_name}_sorted
join -j 1 ./mQTL/GoDMC/${my_list_name}_sorted ./mQTL/GoDMC/snps_36cohorts16_rsid_sorted.txt | weird >./mQTL/GoDMC/snps_36cohorts16_rsid_vs_${my_list_name}
#Output only significant results. A cis pvalue smaller than 1e-8 and a trans pvalue smaller than 1e-14 is deemed as significant in the GoDMC analysis.
cat ./mQTL/GoDMC/snps_36cohorts16_rsid_vs_${my_list_name} | awk -v OFS="\t" '{if ($28=="TRUE" && $10 < 1e-8) print $0  ; else if ($28=="FALSE" && $10 < 1e-14)  print $0}' >./mQTL/GoDMC/snps_36cohorts16_rsid_vs_${my_list_name}_significant
#Use p-values computed using multiplicative random effects meta analysis - more stringent
cat ./mQTL/GoDMC/snps_36cohorts16_rsid_vs_${my_list_name} | awk -v OFS="\t" '{if ($28=="TRUE" && $21 < 1e-8) print $0  ; else if ($28=="FALSE" && $21 < 1e-14)  print $0}' >./mQTL/GoDMC/snps_36cohorts16_rsid_vs_${my_list_name}_significant_mr
#Annotate results with closes TSS and transcript.
module load languages/R-3.5.1-ATLAS-gcc-6.1
Rscript --vanilla $scripts/annotate_cpg.R ./mQTL/GoDMC/snps_36cohorts16_rsid_vs_${my_list_name}_significant_mr
Rscript --vanilla $scripts/annotate_cpg.R ./mQTL/GoDMC/snps_36cohorts16_rsid_vs_${my_list_name}_significant
}

function grep_ishigaki {
my_list_file=$1
my_list_name=$2
#Filter only significant values (Q SNP < 0.05) according to permutation testing.
for my_textfile in ./eQTL/Ishigaki2017/eQTL_exon_level/permutation/B_permutation_q_value_0.5.txt \
./eQTL/Ishigaki2017/eQTL_exon_level/permutation/CD4_permutation_q_value_0.5.txt \
./eQTL/Ishigaki2017/eQTL_exon_level/permutation/CD8_permutation_q_value_0.5.txt \
./eQTL/Ishigaki2017/eQTL_exon_level/permutation/Mono_permutation_q_value_0.5.txt \
./eQTL/Ishigaki2017/eQTL_exon_level/permutation/NK_permutation_q_value_0.5.txt \
./eQTL/Ishigaki2017/eQTL_exon_level/permutation/PB_permutation_q_value_0.5.txt 
do
head -1 $my_textfile >${my_textfile%.txt}_sig.txt
cat $my_textfile | awk -v OFS="\t" '$11 <= 0.05 {print $0}' >>${my_textfile%.txt}_sig.txt
grep_snp_header $my_list_file $my_list_name ${my_textfile%.txt}_sig.txt
done

for my_textfile in ./eQTL/Ishigaki2017/eQTL_gene_level/permutation/B_permutation_q_value_0.5.txt \
./eQTL/Ishigaki2017/eQTL_gene_level/permutation/CD4_permutation_q_value_0.5.txt \
./eQTL/Ishigaki2017/eQTL_gene_level/permutation/CD8_permutation_q_value_0.5.txt \
./eQTL/Ishigaki2017/eQTL_gene_level/permutation/Mono_permutation_q_value_0.5.txt \
./eQTL/Ishigaki2017/eQTL_gene_level/permutation/NK_permutation_q_value_0.5.txt \
./eQTL/Ishigaki2017/eQTL_gene_level/permutation/PB_permutation_q_value_0.5.txt 
do
head -1 $my_textfile >${my_textfile%.txt}_sig.txt
cat $my_textfile | awk -v OFS="\t" '$10 <= 0.05 {print $0}' >>${my_textfile%.txt}_sig.txt
grep_snp_header $my_list_file $my_list_name ${my_textfile%.txt}_sig.txt
done
}

function grep_chosen_files7 {
my_list_file=$1
my_list_name=$2
for my_textfile in ./caQTL/Maurano2015/ng.3432-S5.txt 
do
head -1 $my_textfile >${my_textfile%.txt}_sig.txt
cat $my_textfile | grep "imbalanced_(5%_FDR)" >>${my_textfile%.txt}_sig.txt
grep_snp_header $my_list_file $my_list_name ${my_textfile%.txt}_sig.txt
done

for my_textfile in ./caQTL/Maurano2015/ng.3432-S7.txt \
./promoter-enhancer/Burren2017/13059_2017_1285_MOESM9_ESM.tsv \
./promoter-enhancer/Mumbach2017/ng.3963-S7.*tsv
do 
grep_snp_header $my_list_file $my_list_name $my_textfile 
done

for my_textfile in ./pQTL/Emilsson2018/aaq1327_Excel_tables.9.tsv \
./pQTL/Emilsson2018/aaq1327_Excel_tables.10.tsv \
./pQTL/Emilsson2018/aaq1327_Excel_tables.11.tsv \
./pQTL/Emilsson2018/aaq1327_Excel_tables.12.tsv 
do

grep_snp_header_n $my_list_file $my_list_name $my_textfile 2
done
}

function grep_chosen_files8 {
my_list_file=$1
my_list_name=$2
for my_textfile in ./network/Fagny2017/whole_blood.qtls \
./network/Fagny2017/skin.qtls 
do
head -1 $my_textfile >${my_textfile%.txt}_sig.txt
cat $my_textfile | awk -v OFS="\t" '$10 < 0.05 {print $0}' >>${my_textfile%.txt}_sig.txt
grep_snp_header $my_list_file $my_list_name ${my_textfile%.txt}_sig.txt
done

for my_textfile in ./network/Fagny2017/whole_blood.snps \
./network/Fagny2017/skin.snps
do
grep_snp_header $my_list_file $my_list_name $my_textfile 
done
}


function grep_chosen_files9 {
my_list_file=$1
my_list_name=$2
for my_textfile in ./RNAs/miRNA/Ziebarth2012/SNPs_and_indels_in_miRNA_seeds_human.txt \
./RNAs/miRNA/Ziebarth2012/experimentally_supported_miRSNP_human.txt \
./RNAs/miRNA/Ziebarth2012/target_miRSNP_human_CLASH.txt \
./RNAs/miRNA/Ziebarth2012/Genes_associated_with_human_diseases_traits.txt \
./RNAs/miRNA/Ziebarth2012/target_miRSNP_human.txt
do
grep_snp_header $my_list_file $my_list_name $my_textfile
done
}

function grep_chosen_files10 {
my_list_file=$1
my_list_name=$2
for my_textfile in ./splicing/Xiong2014/hg19_spidex_sig_rsid.txt
do
grep_snp_header $my_list_file $my_list_name $my_textfile
done
for my_textfile in ./rVarBase/dbsnp_info1.txt
do
grep_snp_noheader $my_list_file $my_list_name $my_textfile
done
for my_textfile in ./eQTL/Schmiedel2018/1-s2.0-S009286741831331X-mmc3.1.tsv
do
grep_snp_header_n $my_list_file $my_list_name $my_textfile 7
done
for my_textfile in ./eQTL/Schmiedel2018/*vcf
do
grep_snp_header $my_list_file $my_list_name $my_textfile
done
}

function grep_chosen_files11 {
my_list_file=$1
my_list_name=$2
for my_textfile in ./eQTL/Alasoo2018/41588_2018_46_MOESM3_ESM.5.tsv ./eQTL/Andiappan2015/ncomms8971-s2.1.tsv ./eQTL/Andiappan2015/ncomms8971-s3.1.tsv \
./eQTL/Kim2014/ncomms6236-s4.1.tsv ./eQTL/Kim2014/ncomms6236-s5.1.tsv ./eQTL/Kim2014/ncomms6236-s6.1.tsv ./eQTL/Ye2014/tableS7-cis_eQTL_meta.8.tsv \
./eQTL/Nedelec2016/1-s2.0-S0092867416313071-mmc5.1.tsv ./eQTL/Quach2016/1-s2.0-S009286741631306X-mmc2.1.tsv \
./eQTL/Quach2016/1-s2.0-S009286741631306X-mmc2.3.tsv ./eQTL/Quach2016/1-s2.0-S009286741631306X-mmc3.1.tsv \
./eQTL/Quach2016/1-s2.0-S009286741631306X-mmc4.2.tsv ./eQTL/Quach2016/1-s2.0-S009286741631306X-mmc4.3.tsv \
./eQTL/Lee2014/TableS4_cis_eQTLs_and_reQTLs_from_pooled_analysis.1.tsv \
./eQTL/Lee2014/TableS4_cis_eQTLs_and_reQTLs_from_pooled_analysis.2.tsv \
./eQTL/Lee2014/TableS4_cis_eQTLs_and_reQTLs_from_pooled_analysis.3.tsv \
./eQTL/Lee2014/TableS4_cis_eQTLs_and_reQTLs_from_pooled_analysis.4.tsv \
./eQTL/Lee2014/TableS8_trans_eQTLs_and_reQTLs_from_pooled_analysis.1.tsv \
./eQTL/Lee2014/TableS8_trans_eQTLs_and_reQTLs_from_pooled_analysis.2.tsv \
./eQTL/Lee2014/TableS8_trans_eQTLs_and_reQTLs_from_pooled_analysis.3.tsv \
./eQTL/Lee2014/TableS8_trans_eQTLs_and_reQTLs_from_pooled_analysis.4.tsv \
./eQTL/Pala2017/eQTLs.Fdr0.05.WithConditional.Annot.Recoded.tsv.1k \
./eQTL/Pala2017/isoQTLs.Fdr0.05.WithConditional.Annot.Recoded.tsv.1k \
./eQTL/Naranbhai2015/ncomms8545-s2.1.tsv \
./eQTL/Naranbhai2015/ncomms8545-s4.1.tsv \
./eQTL/Raj2014/tableS4_eu_cd4T_cis_fdr05.tsv \
./eQTL/Raj2014/tableS5_ea_cd4T_cis_fdr05.tsv \
./eQTL/Raj2014/tableS6_aa_cd4T_cis_fdr05.tsv \
./eQTL/Raj2014/tableS7_eu_monocytes_cis_fdr05.tsv \
./eQTL/Raj2014/tableS8_ea_monocytes_cis_fdr05.tsv \
./eQTL/Raj2014/tableS9_aa_monocytes_cis_fdr05.tsv \
./eQTL/Raj2014/tableS11_meta_monocytes_cis_fdr05.tsv \
./eQTL/Raj2014/tableS12_meta_cd4T_cis_fdr05.tsv \
./eQTL/Raj2014/tableS13_trans_monocytes_bf.tsv \
./eQTL/Raj2014/tableS14_trans_cd4T_bf.tsv 
do
grep_snp_header $my_list_file $my_list_name $my_textfile
done
}

#1K r2 > 0.2 grep
for my_f in grep_chosen_files grep_chosen_files2 grep_chosen_files3 grep_chosen_files4 grep_chosen_files5 \
grep_chosen_files6 grep_chosen_files7 grep_chosen_files8 grep_chosen_files8 grep_chosen_files9 grep_chosen_files10 grep_chosen_files11 grep_ishigaki
do
$my_f $data_manipulation/interval_r2_0.2_1k.snps interval_r2_0.2_1k
done
#using the Pvalue column for filtering - in total, 36896 variants have a significant mQTL and they are matched to 6060 CpGs neighbouting 864 genes. An average of 143 mQTLs per variant.
#Using the PValueMR (most stringent), we get 36171 variants have a significant mQTL and they are matched to 5137 CpGs neighbouring 731 genes. An average of 98.5 mQTLs per variant.
grep_godmc $data_manipulation/interval_r2_0.2_1k.godmc interval_r2_0.2_1k
#Include only significant p-values for Suhre et al.
awk -v OFS="\t" '$NF < 10e-11 {print $0}' ./pQTL/Suhre2017/all_suhre2017.txt_vs_interval_r2_0.2_1k > ./pQTL/Suhre2017/all_suhre2017.txt_vs_interval_r2_0.2_1k.sig
#Merge two tables from Fagny et al. 2017
cd $datasets/network/Fagny2017
qsub $scripts/sub_fagny2017_merge.sh
cd -

#General overview of the number of hits per SNP.
for my_file in ./*/*/*_vs_interval_r2_0.2_1k*
do
file_basename=$(basename $my_file)
python $scripts/annotation_overview_counts.py $data_manipulation/interval_r2_0.2_1k.snps $my_file ${file_basename}.counts
done

#Combine all the individuals tables together.
Rscript --vanilla $scripts/merge_counts.R interval_r2_0.2.summed

#Create a version of the file above, where only one hit is counter per analysed file.
awk -v OFS="\t" '{for(i=1;i<=NF;i++) if($i>1) {$i="1"}; print}' interval_r2_0.2.summed >interval_r2_0.2.summed.single