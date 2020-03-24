#!/bin/bash
HOME=/panfs/panasas01/sscm/qh18484
analysis=$HOME/analysis/colocalization/coloc/eqtlgen/euro
scripts=$HOME/bin/eczema_gwas_fu/colocalization/coloc
gwas=/panfs/panasas01/sscm/qh18484/data/gwas/paternoster2015
eqtlgen=/panfs/panasas01/sscm/qh18484/data/eqtl/eQTLgen
utils=$HOME/bin/eczema_gwas_fu/utils

eqtl_f=eQTLsFDR-ProbeLevel.txt.gz
gwas_f=results.euro.pval.eqtlgen.ea

cd $analysis

#Generate Euro GWAS file with rsids matching those from the eQTLgen study.
python $utils/update_rsid.py --bim $HOME/analysis/colocalization/coloc/eqtlgen/published/diagnostics/${eqtl_f%.txt.gz}_unique.bim \
--tab $gwas/results.euro.pval.tsv --head Y --chrom 2 --pos 3 --ref 4 --alt 5 >$analysis/results.euro.pval.eqtlgen

#Harmonize the effect sizes so that they are relevant to the same allele in both files.
python $utils/harmonize_beta.py --tab $analysis/results.euro.pval.eqtlgen --ref $eqtlgen/$eqtl_f --header_tab Y --header_ref Y \
--rsid_tab 12 --rsid_ref 2 --effect_tab 4 --alt_tab 5 --effect_ref 10 --beta_tab 8 --zscore_tab 10 --out $analysis/$gwas_f

#Export eQTLgen results in a given interval.
for my_int in 10 100 250 1000
do
python $scripts/generate_eqtl_input_pvals_zscore.py --tab $eqtlgen/$eqtl_f \
--proces $gwas/paternoster_2015_index_snps_sorted.txt --interval $my_int --chrom 3 --pos 4 --pval 1 --ident 2 --zscore 11
done

#Export GWAS results in a given interval.
for my_int in 10 100 250 1000
do
python $scripts/generate_gwas_input.py --tab $analysis/$gwas_f \
--proces $gwas/paternoster_2015_index_snps_sorted.txt --interval $my_int --chrom 2 --pos 3 --eaf 7 --ident 12 
done

#Run coloc on all input files. Use two analyses types: just with p-values, and with betas.
#Noteby only top eQTL for each SNP is considered so our gene of interest may be hidden down the list?
for my_input in *.eqtl
do
Rscript $scripts/coloc_eqtlgen2.R ${my_input%.eqtl}.gwas $my_input 
done

#Subset to only results with robust support for colocalisation (p > 0.75)
for my_input in *.colocp
do
awk '{if (NR==1) ; else if ($15 > 0.75 && $15 <= 1) print $0}' $my_input >${my_input%.coloc}.sig
done

#Subset to only results with significant support for colocalisation (p > 0.95 )
for my_input in *.colocp
do
awk -F '\t' '{if (NR==1) ; else if ($15 > 0.95 && $15 <= 1) print $0}' $my_input >${my_input%.coloc}.sig95
done

#Annotate all the loci with significant probability of colocalisation (p> 0.75) between SNP in the GWAS interval and eQTL 
#with the transcript ID potentially affected
for a in *.sig
do
python $scripts/annotate_eqtlgen.py $a $eqtlgen/eQTLsFDR-ProbeLevel.txt.gz_10e-20.annotable
done

#Same, but for p>0.95
for a in *.sig95
do
python $scripts/annotate_eqtlgen.py $a $eqtlgen/eQTLsFDR-ProbeLevel.txt.gz_10e-20.annotable
done

#Annotate loci with Ensembl.
qsub $scripts/sub_ensembl_coloc_sig.sh

#Plot all results with LocusZoom
for my_snps in rs61813875 rs10791824 rs6419573 rs12188917 rs2212434 rs2918307 rs4809219 rs4713555 rs2041733 rs2944542 rs145809981 rs6827756 rs12730935 rs4312054 rs1249910 rs2592555 rs7127307 rs2038255 rs6602364 rs6473227 rs7512552 rs112111458 rs7625909 rs4643526 rs10199605 rs10214237 rs12951971 rs1057258 rs6872156 rs7016497 rs2905493 rs1799986 rs2227483 rs7146581 rs11657987 rs77714197 rs3917265 rs13152362 rs4705962 rs2218565 rs7512552 rs12730935 rs12295535 rs11923593 rs2143950 rs17881320
do
for my_int in 10 100 250 1000
do
a=${my_snps}_${my_int}kbp.colocp
qsub -v input_file=$analysis/$a,snp_id="snp",pval="SNP.PP.H4",my_ref=${my_snps},my_flank=${my_int}kb,my_prefix=${my_snps}_${my_int}kbp.colocp $utils/sub_locus_zoom.sh 
done
done
#Same, but for p>0.95

#Convert to PNG.
convert -verbose -density 500 "${my_pdf}" "${my_pdf%.*}.png"
#Concatenate all annotated significant results
for a in *sig.ensembl.gene
do
echo ${a%.colocp.sig95.ensembl.gene} >>eqltgen.euro.summary.sig
cat $a >> eqltgen.euro.summary.sig
done

#Grep for nearest genes to the index SNP.
while read line
do
set -- $line
gene=$1
grep $gene eqltgen.euro.summary.sig
done < $gwas/paternoster2015_gene_list

#Convert chosen PDF LocusZoom plots to PNG.
for my_pdf in rs10214237_10kbp.colocp_rs10214237.pdf rs17881320_10kbp.colocp_rs17881320.pdf
do
convert -verbose -density 500 "${my_pdf}" "${my_pdf%.*}.png" 
done 

#Identify the genes with overall significant PP.H4 calculated using p-vals.
for a in *totalp
do
awk 'NR == 7 && $2 > 0.65 {print FILENAME}' $a
done
#11 genes total.
ls rs12295535_1000kbp.totalp
ls rs12295535_100kbp.totalp
ls rs12295535_10kbp.totalp
ls rs12295535_250kbp.totalp
ls rs2038255_100kbp.totalp
ls rs2038255_10kbp.totalp
ls rs2143950_100kbp.totalp
ls rs2143950_10kbp.totalp
ls rs2212434_10kbp.totalp
ls rs2218565_100kbp.totalp
ls rs2218565_10kbp.totalp
ls rs2218565_250kbp.totalp
ls rs2592555_1000kbp.totalp
ls rs2592555_100kbp.totalp
ls rs2592555_10kbp.totalp
ls rs2592555_250kbp.totalp
ls rs4705962_10kbp.totalp
ls rs4809219_1000kbp.totalp
ls rs4809219_100kbp.totalp
ls rs4809219_10kbp.totalp
ls rs4809219_250kbp.totalp
ls rs6419573_100kbp.totalp
ls rs6419573_10kbp.totalp
ls rs6872156_10kbp.totalp
ls rs7512552_10kbp.totalp

for a in *totalb
do
awk 'NR == 7 && $2 > 0.65 {print FILENAME}' $a
done

#Identify the genes with overall significant PP.H4 calculated using betas and variance of beta.
#10 genes total.
ls rs12295535_1000kbp.totalb
ls rs12295535_100kbp.totalb
ls rs12295535_10kbp.totalb
ls rs12295535_250kbp.totalb
ls rs2038255_1000kbp.totalb
ls rs2038255_100kbp.totalb
ls rs2038255_10kbp.totalb
ls rs2038255_250kbp.totalb
ls rs2143950_1000kbp.totalb
ls rs2143950_100kbp.totalb
ls rs2143950_10kbp.totalb
ls rs2143950_250kbp.totalb
ls rs2212434_10kbp.totalb
ls rs2218565_100kbp.totalb
ls rs2218565_10kbp.totalb
ls rs2218565_250kbp.totalb
ls rs2592555_1000kbp.totalb
ls rs2592555_100kbp.totalb
ls rs2592555_10kbp.totalb
ls rs2592555_250kbp.totalb
ls rs4705962_10kbp.totalb
ls rs4809219_1000kbp.totalb
ls rs4809219_100kbp.totalb
ls rs4809219_10kbp.totalb
ls rs4809219_250kbp.totalb
ls rs6872156_10kbp.totalb
ls rs7512552_10kbp.totalb
