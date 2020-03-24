#!/bin/bash
HOME=/panfs/panasas01/sscm/qh18484
analysis=$HOME/analysis/annotation/giggle
scripts=$HOME/bin/eczema_gwas_fu/annotation/giggle
giggle=$HOME/bin/giggle
data_manipulation=$HOME/analysis/annotation/data_manipulation
datasets=$HOME/working/data/Datasets

cd $analysis

#!! remember that giggle requires tab-delimited fields, space will not do

#Gzip our intervals of interest.
bgzip -c $data_manipulation/interval_r2_0.2_1k.bed > interval_r2_0.2_1k.bed.gz
bgzip -c $data_manipulation/interval_r2_0.2_1k_individual.bed > interval_r2_0.2_1k_individual.bed.gz
bgzip -c $data_manipulation/paternoster_2015_index_snps_sorted_1Mbp.bed > paternoster_2015_index_snps_sorted_1Mbp.bed.gz
bgzip -c $data_manipulation/paternoster_2015_index_snps_sorted_3Mbp.bed > paternoster_2015_index_snps_sorted_3Mbp.bed.gz

$giggle/bin/giggle index -i "$giggle/roadmap_sort/*.gz" -o $giggle/roadmap_sort_index -s 
$giggle/bin/giggle index -i "$giggle/fantom_sort/*.gz" -o $giggle/fantom_sort_index -s 

for a in interval_r2_0.2_1k.bed.gz paternoster_2015_index_snps_sorted_1Mbp.bed.gz paternoster_2015_index_snps_sorted_3Mbp.bed.gz
do
$giggle/bin/giggle search -i $giggle/roadmap_sort_index -q $a -s > ${a%.bed.gz}_roadmap.overview
$giggle/bin/giggle search -i $giggle/fantom_sort_index -q $a -s > ${a%.bed.gz}_fantom.overview
done 

for a in interval_r2_0.2_1k.bed.gz paternoster_2015_index_snps_sorted_1Mbp.bed.gz paternoster_2015_index_snps_sorted_3Mbp.bed.gz
do
$giggle/bin/giggle search -i $giggle/roadmap_sort_index -q $a -o > ${a%.bed.gz}_roadmap.o
$giggle/bin/giggle search -i $giggle/fantom_sort_index -q $a -o > ${a%.bed.gz}_fantom.o
done 

for a in interval_r2_0.2_1k.bed.gz paternoster_2015_index_snps_sorted_1Mbp.bed.gz paternoster_2015_index_snps_sorted_3Mbp.bed.gz
do
$giggle/bin/giggle search -i $giggle/roadmap_sort_index -q $a -v > ${a%.bed.gz}_roadmap.v
$giggle/bin/giggle search -i $giggle/fantom_sort_index -q $a -v > ${a%.bed.gz}_fantom.v
done 

for a in interval_r2_0.2_1k.bed.gz paternoster_2015_index_snps_sorted_1Mbp.bed.gz paternoster_2015_index_snps_sorted_3Mbp.bed.gz
do
$giggle/bin/giggle search -i $giggle/roadmap_sort_index -q $a -v -o > ${a%.bed.gz}_roadmap.ov
$giggle/bin/giggle search -i $giggle/fantom_sort_index -q $a -v -o > ${a%.bed.gz}_fantom.ov
done 

#Prepare results from Jarre et al (2016) into BED format and looks for overlap
#bait
cut -f1,2,3,9 $datasets/promoter-enhancer/Javierre2016/ActivePromoterEnhancerLinks.tsv | tail -n +2 >ActivePromoterEnhancerLinks.bait
sort ActivePromoterEnhancerLinks.bait | uniq >temp
mv temp ActivePromoterEnhancerLinks.bait.bed
#oe
cut -f5,6,7,9 $datasets/promoter-enhancer/Javierre2016/ActivePromoterEnhancerLinks.tsv | tail -n +2 >ActivePromoterEnhancerLinks.oe
sort ActivePromoterEnhancerLinks.oe | uniq >temp
mv temp ActivePromoterEnhancerLinks.oe.bed

cut -f1,2,3 $datasets/promoter-enhancer/Javierre2016/PCHiC_peak_matrix_cutoff5.tsv | tail -n +2 >PCHiC_peak_matrix_cutoff5.bait
sort PCHiC_peak_matrix_cutoff5.bait | sed 's/^/chr/' | uniq >temp
mv temp PCHiC_peak_matrix_cutoff5.bait.bed

cut -f6,7,8 $datasets/promoter-enhancer/Javierre2016/PCHiC_peak_matrix_cutoff5.tsv | tail -n +2 >PCHiC_peak_matrix_cutoff5.oe
sort PCHiC_peak_matrix_cutoff5.oe | sed 's/^/chr/' | uniq >temp
mv temp PCHiC_peak_matrix_cutoff5.oe.bed 

cut -f1,2,3,4,5,9- $datasets/promoter-enhancer/Javierre2016/PCHiC_vs_rCHiC_peak_matrix.tsv | tail -n +2 >PCHiC_vs_rCHiC_peak_matrix.bait
sort PCHiC_vs_rCHiC_peak_matrix.bait | sed 's/^/chr/' | uniq >temp
mv temp PCHiC_vs_rCHiC_peak_matrix.bait.bed

cat $datasets/promoter-enhancer/Javierre2016/PCHiC_vs_rCHiC_peak_matrix.tsv | awk -v OFS="\t" '{print $6, $7, $8, $4, $5, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19}'| tail -n +2 >PCHiC_vs_rCHiC_peak_matrix.oe
sort PCHiC_vs_rCHiC_peak_matrix.oe | sed 's/^/chr/' | uniq >temp
mv temp PCHiC_vs_rCHiC_peak_matrix.oe.bed

#Mumbach et al. (2017)
for a in $datasets/promoter-enhancer/Mumbach2017/ng.3963-S4.*tsv
do
cut -f 1,2,3,7- $a | tail -n +2 | sed 's/^/chr/' | sort | uniq >${a%.tsv}_a.bed
cut -f 4- $a | tail -n +2 | sed 's/^/chr/' | sort | uniq  >${a%.tsv}_b.bed
done

for a in $datasets/promoter-enhancer/Mumbach2017/ng.3963-S5.*tsv $datasets/promoter-enhancer/Mumbach2017/ng.3963-S6.*tsv $datasets/promoter-enhancer/Mumbach2017/ng.3963-S7.*tsv
do
awk -v OFS="\t" '$NF < 0.05 {print $0}' $a | cut -f 1,2,3,4,8- | tail -n +2 | sort | uniq >${a%.tsv}_x.bed
awk -v OFS="\t" '$NF < 0.05 {print $0}' $a | cut -f 5- | tail -n +2 | sort | uniq >${a%.tsv}_y.bed
done

cp -r $datasets/promoter-enhancer/Mumbach2017/*bed ./

cat $datasets/promoter-enhancer/Rinaldi2016/1-s2.0-S1934590916301953-mmc2.1.tsv | awk -v OFS="\t" '$NF == "yes" {print $14}' | sed 's/:/\t/' | sed 's/-/\t/' >1-s2.0-S1934590916301953-mmc2.1.bed

cat $datasets/promoter-enhancer/Rinaldi2016/1-s2.0-S1934590916301953-mmc3.2.tsv | tail -n +5 | cut -f1,2,3 >1-s2.0-S1934590916301953-mmc3.2_EpiSC_DNMT3A_peaks.bed 

cat $datasets/promoter-enhancer/Rinaldi2016/1-s2.0-S1934590916301953-mmc3.2.tsv | tail -n +5 | cut -f5,6,7 >1-s2.0-S1934590916301953-mmc3.2_EpiSC_DNMT3B_peaks.bed 

cat $datasets/promoter-enhancer/Rinaldi2016/1-s2.0-S1934590916301953-mmc3.2.tsv | tail -n +5 | cut -f10,11,12 >1-s2.0-S1934590916301953-mmc3.2_Diff_DNMT3A_peaks.bed 

cat $datasets/promoter-enhancer/Rinaldi2016/1-s2.0-S1934590916301953-mmc3.2.tsv | tail -n +5 | cut -f15,16,17 >1-s2.0-S1934590916301953-mmc3.2_Diff_DNMT3B_peaks.bed 

cat $datasets/promoter-enhancer/Rinaldi2016/1-s2.0-S1934590916301953-mmc4.1.tsv | tail -n +5 | cut -f1,2,3 >1-s2.0-S1934590916301953-mmc4.1_EpiSC_H3K4me1.bed 
cat $datasets/promoter-enhancer/Rinaldi2016/1-s2.0-S1934590916301953-mmc4.1.tsv | tail -n +5 | cut -f6,7,8 >1-s2.0-S1934590916301953-mmc4.1_EpiSC_H3K4m3.bed 
cat $datasets/promoter-enhancer/Rinaldi2016/1-s2.0-S1934590916301953-mmc4.1.tsv | tail -n +5 | cut -f11,12,13 >1-s2.0-S1934590916301953-mmc4.1_EpiSC_H3K27ac.bed 
cat $datasets/promoter-enhancer/Rinaldi2016/1-s2.0-S1934590916301953-mmc4.1.tsv | tail -n +5 | cut -f16,17,18 >1-s2.0-S1934590916301953-mmc4.1_Diff_H3K4m1.bed 
cat $datasets/promoter-enhancer/Rinaldi2016/1-s2.0-S1934590916301953-mmc4.1.tsv | tail -n +5 | cut -f21,22,23 >1-s2.0-S1934590916301953-mmc4.1_Diff_H3K4m3.bed 
cat $datasets/promoter-enhancer/Rinaldi2016/1-s2.0-S1934590916301953-mmc4.1.tsv | tail -n +5 | cut -f26,27,28 >1-s2.0-S1934590916301953-mmc4.1_Diff_H3K27ac.bed 

cat $datasets/promoter-enhancer/Rinaldi2016/1-s2.0-S1934590916301953-mmc4.2.tsv | tail -n +6 | cut -f1,2,3 >1-s2.0-S1934590916301953-mmc4.2.a.bed
cat $datasets/promoter-enhancer/Rinaldi2016/1-s2.0-S1934590916301953-mmc4.2.tsv | tail -n +6 | cut -f6,7,8 >1-s2.0-S1934590916301953-mmc4.2.b.bed
cat $datasets/promoter-enhancer/Rinaldi2016/1-s2.0-S1934590916301953-mmc4.2.tsv | tail -n +6 | cut -f14,15,16 >1-s2.0-S1934590916301953-mmc4.2.c.bed
cat $datasets/promoter-enhancer/Rinaldi2016/1-s2.0-S1934590916301953-mmc4.2.tsv | tail -n +6 | cut -f20,21,22 >1-s2.0-S1934590916301953-mmc4.2.d.bed
cat $datasets/promoter-enhancer/Rinaldi2016/1-s2.0-S1934590916301953-mmc4.2.tsv | tail -n +6 | cut -f26,27,28 >1-s2.0-S1934590916301953-mmc4.2.e.bed

cat $datasets/promoter-enhancer/Rinaldi2016/1-s2.0-S1934590916301953-mmc4.3.tsv | tail -n +6 | cut -f19,20,21 >1-s2.0-S1934590916301953-mmc4.3.a.bed
cat $datasets/promoter-enhancer/Rinaldi2016/1-s2.0-S1934590916301953-mmc4.3.tsv | tail -n +6 | cut -f31,32,33 >1-s2.0-S1934590916301953-mmc4.3.b.bed

cat $datasets/promoter-enhancer/Rinaldi2016/stem_2032_mmc5.1.tsv | tail -n +10 | cut -f2,3,4,5,8 >stem_2032_mmc5.1.a.bed
cat $datasets/promoter-enhancer/Rinaldi2016/stem_2032_mmc5.1.tsv | tail -n +10 | cut -f15,16,17,18,21 >stem_2032_mmc5.1.b.bed

cat $datasets/promoter-enhancer/Rinaldi2016/stem_2032_mmc5.2.tsv | tail -n +9 | cut -f2,3,4 >stem_2032_mmc5.2.a.bed
cat $datasets/promoter-enhancer/Rinaldi2016/stem_2032_mmc5.2.tsv | tail -n +9 | cut -f14,15,16 >stem_2032_mmc5.2.b.bed

cat $datasets/promoter-enhancer/Burren2017/13059_2017_1285_MOESM3_ESM | tail -n +2 | awk -v OFS="\t" '$8=="oe" && $15 < 0.05 {print $0}' | weird | sort | uniq >13059_2017_1285_MOESM3_ESM.oe.bed
cat $datasets/promoter-enhancer/Burren2017/13059_2017_1285_MOESM3_ESM | tail -n +2 | awk -v OFS="\t" '$8=="bait" && $15 < 0.05 {print $0}'  | weird | sort | uniq >13059_2017_1285_MOESM3_ESM.bait.bed
cat $datasets/promoter-enhancer/Burren2017/13059_2017_1285_MOESM5_ESM | tail -n +2  | awk -v OFS="\t" '{print $4, $5, $5+$6}' | sort | uniq >13059_2017_1285_MOESM5_ESM.bait.bed
cat $datasets/promoter-enhancer/Burren2017/13059_2017_1285_MOESM5_ESM | tail -n +2  | awk -v OFS="\t" '{print $7, $8, $9+$8}' | sort | uniq >13059_2017_1285_MOESM5_ESM.oe.bed
cat $datasets/promoter-enhancer/Javierre2016/RegBuild_BLUEPRINT_annotations.txt | tail -n +2  | cut -f2-4,9- | sort | uniq >RegBuild_BLUEPRINT_annotations.bait.bed
cat $datasets/promoter-enhancer/Javierre2016/RegBuild_BLUEPRINT_annotations.txt | tail -n +2  | cut -f6- | sort | uniq >RegBuild_BLUEPRINT_annotations.other.bed
cat $datasets/promoter-enhancer/Mifsud2015/TS5_GM12878_promoter-other_significant_interactions.txt | tail -n +2 | cut -f 1-6 >TS5_GM12878_promoter-other_significant_interactions.x.bed
cat $datasets/promoter-enhancer/Mifsud2015/TS5_GM12878_promoter-other_significant_interactions.txt | tail -n +2 | awk -v OFS="\t" '{print $7, $8, $9, $4, $5, $6}' >TS5_GM12878_promoter-other_significant_interactions.y.bed
cat $datasets/promoter-enhancer/Mifsud2015/TS5_CD34_promoter-other_significant_interactions.txt | tail -n +2 | cut -f 1-6 >TS5_CD34_promoter-other_significant_interactions.x.bed
cat $datasets/promoter-enhancer/Mifsud2015/TS5_CD34_promoter-other_significant_interactions.txt | tail -n +2 | awk -v OFS="\t" '{print $7, $8, $9, $4, $5, $6}'>TS5_CD34_promoter-other_significant_interactions.y.bed
cat $datasets/promoter-enhancer/Mifsud2015/TS5_GM12878_promoter-promoter_significant_interactions.txt | tail -n +2 | cut -f 1-6,13- >TS5_GM12878_promoter-promoter_significant_interactions.x.bed
cat $datasets/promoter-enhancer/Mifsud2015/TS5_GM12878_promoter-promoter_significant_interactions.txt | tail -n +2 | cut -f 7- >TS5_GM12878_promoter-promoter_significant_interactions.y.bed
cat $datasets/promoter-enhancer/Mifsud2015/TS5_CD34_promoter-promoter_significant_interactions.txt | tail -n +2 | cut -f 1-6,13- >TS5_CD34_promoter-promoter_significant_interactions.x.bed
cat $datasets/promoter-enhancer/Mifsud2015/TS5_CD34_promoter-promoter_significant_interactions.txt | tail -n +2 | cut -f 7- >TS5_CD34_promoter-promoter_significant_interactions.y.bed
for a in interval_r2_0.2_1k_individual.bed.gz_vs_TS5_*_promoter-promoter_significant_interactions.*.bed_index_processed.bed
do
sed -i 's/|/;/g' $a
done

for a in interval_r2_0.2_1k_individual.bed.gz_vs_TS5_*_promoter-other_significant_interactions.*.bed_index_processed.bed
do
sed -i 's/|/;/g' $a
done

cat $datasets/promoter-enhancer/Gasperini2019/mmc2.2-3_combined.txt | tail -n +2  | awk -v OFS="\t" '{print $9, $10, $11, $1, $2, $3, $4, $5, $6, $7, $8}' >mmc2.2-3_combined.bed
cat $datasets/promoter-enhancer/Rubin2017/ng.3935-S5.1.tsv | tail -n +4 | cut -f1-3 >ng.3935-S5.1.x.bed
cat $datasets/promoter-enhancer/Rubin2017/ng.3935-S5.1.tsv | tail -n +4 | cut -f4-6 >ng.3935-S5.1.y.bed
cat $datasets/promoter-enhancer/Rubin2017/ng.3935-S5.2.tsv | tail -n +4 | cut -f1-3 >ng.3935-S5.2.x.bed
cat $datasets/promoter-enhancer/Rubin2017/ng.3935-S5.2.tsv | tail -n +4 | cut -f4-6 >ng.3935-S5.2.y.bed
cat $datasets/promoter-enhancer/Rubin2017/ng.3935-S5.3.tsv | tail -n +4 | cut -f1-3 >ng.3935-S5.3.x.bed
cat $datasets/promoter-enhancer/Rubin2017/ng.3935-S5.3.tsv | tail -n +4 | cut -f4-6 >ng.3935-S5.3.y.bed
cat $datasets/promoter-enhancer/Rubin2017/ng.3935-S5.4.tsv | tail -n +4 | cut -f1-3 >ng.3935-S5.4.x.bed
cat $datasets/promoter-enhancer/Rubin2017/ng.3935-S5.4.tsv | tail -n +4 | cut -f4-6 >ng.3935-S5.4.y.bed

cat $datasets/promoter-enhancer/Freire-Pritchett2017/elife-21926-supp1-v2.txt | tail -n +2 | cut -f 1-3,9- >elife-21926-supp1-v2.bait.bed
cat $datasets/promoter-enhancer/Freire-Pritchett2017/elife-21926-supp1-v2.txt | tail -n +2 | cut -f 5-7,9- >elife-21926-supp1-v2.pir.bed
cat $datasets/promoter-enhancer/Freire-Pritchett2017/elife-21926-supp2-v2.txt | tail -n +2 | sed 's/[:-]/\t/g' >elife-21926-supp2-v2.bed
cp $datasets/TADs/Freire-Pritchett2017/hESC_TADs_delta2.0.bed ./

cat $datasets/promoter-enhancer/ReMap2018/1k_r2_02_min10percent_peaks/intersectBed_20181109120309000000_2060987663.bed | cut -f 1-4 | sort | uniq >ReMap2018.bed

cp $datasets/promoter-enhancer/Wang2015/pnas.1507253112.sd01_hg19.bed ./
cp $datasets/caQTL/Lander2017/*.bed ./
cp -r $datasets/RNAs/lncRNA/Chen2012/human_hg37.bed ./
cp -r $datasets/RNAs/lncRNA/Xie2014/NONCODEv5_hg19.lncAndGene.bed ./
cp -r $datasets/RNAs/piRNA/Wang2018/piR_hg19_sort.bed ./
cp -r $datasets/RNAs/miRNA/Rie2017/nbt.3947-S8.txt nbt.3947-S8.bed

cut -f5,6,7,11 $datasets/promoter-enhancer/Ziebarth2013/CTCFBSDB_all_exp_sites_Sept12_2012_human_hg19.txt >CTCFBSDB_all_exp_sites_Sept12_2012_human_hg19.bed

cut -f5,6,7 $datasets/promoter-enhancer/Ziebarth2013/CTCFBSDB_all_exp_sites_Sept12_2012_human_hg18.txt >CTCFBSDB_all_exp_sites_Sept12_2012_human_hg18.bed
tail -n +2 $datasets/promoter-enhancer/Wang2017/csre.tab >csre.bed
cp -r $datasets/TADs/Dixon2012/hESC_total.combined.domain hESC_total.combined.domain.bed
cp -r $datasets/TADs/Dixon2012/IMR90_total.combined.domain IMR90_total.combined.domain.bed
cp -r $datasets/TADs/ENCODE/*bed ./
cp -r $datasets/TADs/Harmston2017/41467_2017_524_MOESM2_ESM.bed ./
cp -r $datasets/TADs/Rao2014/GM12878_CTCF_orientation.bed ./
cp -r $datasets/hQTL/Pelikan2018/pelikan2018_supp_data3.bed ./
cp -r $datasets/TADs/Javierre2016/*bed ./
cp -r $datasets/promoter-enhancer/Mifsud2015/mifsud2015_supp_tab2_LCL_TADs.bed ./
cp -r $datasets/promoter-enhancer/Wang2018/41467_2018_7746_MOESM6_ESM.txt ./41467_2018_7746_MOESM6_ESM.bed
cp -r $datasets/promoter-enhancer/Wang2018/41467_2018_7746_MOESM4_ESM.txt ./41467_2018_7746_MOESM4_ESM.bed

tail -n +2 $datasets/promoter-enhancer/Fulco2019/media-3.1.tsv >media-3.1.bed

#Create index and run Giggle
for a in *.bed
do 
qsub -v input=$a,output=${a}_index $scripts/sub_giggle_index.sh 
done

for a in *.bed
do
qsub -v input=${a}_index,query=interval_r2_0.2_1k.bed.gz $scripts/sub_giggle.sh
qsub -v input=${a}_index,query=interval_r2_0.2_1k_individual.bed.gz $scripts/sub_giggle.sh
done

#Turn into per-SNP analysis.
for a in interval_r2_0.2_1k_individual.bed.gz_vs_*.gz
do 
python $scripts/process_giggle_output.py $a ${a%.gz}_processed.bed
done

#Summarise number of matches per SNP per file.
for my_file in *_processed.bed
do
python $HOME/bin/eczema_gwas_fu/annotation/QTL-overlap/annotation_overview_counts.py \
$data_manipulation/interval_r2_0.2_1k.snps $my_file ${my_file}.counts
done

#Merge count tables.
Rscript --vanilla $HOME/bin/eczema_gwas_fu/annotation/QTL-overlap/merge_counts.R interval_r2_0.2.summed

#Create a version of the file above, where only one hit is counter per analysed file.
awk -v OFS="\t" '{for(i=1;i<=NF;i++) if($i>1) {$i="1"}; print}' interval_r2_0.2.summed >interval_r2_0.2.summed.single