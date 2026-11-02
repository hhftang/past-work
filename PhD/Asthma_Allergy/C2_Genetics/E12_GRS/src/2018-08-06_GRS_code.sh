#!/bin/bash

# ==============================================================================
# GET GRS SCORES      
# ==============================================================================

# rawsummary statistics in 
# /projects/GenomicPrediction/GRS_resources/raw_SummaryStatistics

INDIR=/projects/GenomicPrediction/GRS_resources/raw_SummaryStatistics/
OUTDIR=/home/htang/projects/htang/CAS/C2_genetics/GRS
BFILE=/projects/qinqinhuang/CAS/Genotype_Data/process_data/Michigan_imputation_server_results/cas_auto_sub_filt-updated_auto_chr_2018_01_11_01_15/filtered_geno/cas_michigan_imp_filtered

# /projects/Microbiome/CAS_GRS/ for completed CAS scores from Shumei 
# for selected traits

# ------------------------------------------------------------------------------
# ASTHMA
# ------------------------------------------------------------------------------

# Format summary statistics file 
# RSID, alternate allele, European_ancestry_beta_fixed, beta_random, Chr:position, Reference allele
cat "$INDIR"/Asthma/TAGC_Multiancestry_and_European-Ancestry_Meta-analyses_Results.tsv |\
	awk '{ printf $2"\t"$5"\t"$15"\t"$18"\t"$1":"$3"\t"$4"\n"}' > \
	"$OUTDIR"/sumstats/TAGC_European.tsv

# Take the subset of SNPs that matches UKB r2<0.9 LD-thinned set
grep -w -F -f \
	"$OUTDIR"/ukb_imp_allchr_v2_maf0.001_info0.01_whitebritish_tuning_r2_0.9_noambig.prune.in \
	"$OUTDIR"/sumstats/TAGC_European.tsv > \
	"$OUTDIR"/sumstats/asthma_grs_weights.txt

# Check for duplicate rsID
cat "$OUTDIR"/sumstats/asthma_grs_weights.txt | awk '{print $1}'| sort | uniq -d > \
	"$OUTDIR"/sumstats/asthma_duplicate_SNPS.txt

# Calculate GRS profile for CAS based on best SNP list
# Arguments for score represent: file, column of variant ID, allele ID, beta/score 
plink --bfile $BFILE --score "$OUTDIR"/sumstats/asthma_grs_weights.txt 5 2 3 sum \
 --memory 64000 --out "$OUTDIR"/asthma

# ------------------------------------------------------------------------------
# Allergic diseases
# ------------------------------------------------------------------------------

# Format summary statistics
# RSID, EFFECT_ALLELE, BETA SNP (Chr:position) OTHER_ALLELE
cat "$INDIR"/Allergic_disease/SHARE-without23andMe.LDSCORE-GC.SE-META.v0 | \
	awk '{ printf $10"\t"toupper($4)"\t"$6"\t"$1"\t"$5"\n"}' > \
	"$OUTDIR"/sumstats/Allergic_disease.tsv

# Take the subset of SNPs that matches UKB r2<0.9 LD-thinned set
grep -w -F -f \
	"$OUTDIR"/ukb_imp_allchr_v2_maf0.001_info0.01_whitebritish_tuning_r2_0.9_noambig.prune.in \
	"$OUTDIR"/sumstats/Allergic_disease.tsv > \
	"$OUTDIR"/sumstats/allergic_grs_weights.txt

# Check for duplicate rsID
cat "$OUTDIR"/sumstats/allergic_grs_weights.txt | awk '{print $1}'| sort | uniq -d > \
	"$OUTDIR"/sumstats/allergic_duplicate_SNPS.txt

# Rename one SNP
sed -i -e 's/esv3630865;//' "$OUTDIR"/sumstats/allergic_grs_weights.txt

#Calculate GRS profile for CAS based on best SNP list

plink --bfile $BFILE --score "$OUTDIR"/sumstats/allergic_grs_weights.txt 4 2 3 sum \
	--memory 64000 --out "$OUTDIR"/allergic

# ------------------------------------------------------------------------------
# Allergic rhinitis
# ------------------------------------------------------------------------------

# Check for duplicate rsID
cat "$OUTDIR"/sumstats/ar_grs_weights_V2.txt | awk '{print $1}'| sort | uniq -d > \
	"$OUTDIR"/sumstats/ar_duplicate_SNPS.txt

#Calculate GRS profile for CAS based on best SNP list

plink --bfile $BFILE --score "$OUTDIR"/sumstats/ar_grs_weights_V2.txt 4 2 3 sum \
	--memory 64000 --out "$OUTDIR"/ar

# ------------------------------------------------------------------------------
# COPD
# ------------------------------------------------------------------------------

# Format summary statistics
# RSID, EFFECT_ALLELE, BETA SNP (Chr:position) OTHER_ALLELE

tail -n +2 "20002_1112.assoc.tsv" | \
	awk 'BEGIN {FS = ":"; OFS = "\t"} {print $1, $2, $3, $4;}' | \
	awk 'BEGIN {FS="\t"} { printf $5"\t"$4"\t"$9"\t"$1":"$2"\t"$3"\n";}' > \
	"$OUTDIR"/sumstats/COPD.tsv

head -n 1 "$OUTDIR"/sumstats/Allergic_disease.tsv | cat - "$OUTDIR"/sumstats/COPD.tsv > \
	temp && mv temp "$OUTDIR"/sumstats/COPD.tsv 

# Take the subset of SNPs that matches UKB r2<0.9 LD-thinned set
grep -w -F -f \
	"$OUTDIR"/ukb_imp_allchr_v2_maf0.001_info0.01_whitebritish_tuning_r2_0.9_noambig.prune.in \
	"$OUTDIR"/sumstats/COPD.tsv > \
	"$OUTDIR"/sumstats/COPD_grs_weights.txt

# Check for duplicate rsID
cat "$OUTDIR"/sumstats/COPD_grs_weights.txt | awk '{print $1}'| sort | uniq -d > \
	"$OUTDIR"/sumstats/COPD_duplicate_SNPS.txt

# Rename one SNP
sed -i -e 's/esv3630865;//' "$OUTDIR"/sumstats/COPD_grs_weights.txt

#Calculate GRS profile for CAS based on best SNP list

plink --bfile $BFILE --score "$OUTDIR"/sumstats/COPD_grs_weights.txt 4 2 3 sum \
	--memory 64000 --out "$OUTDIR"/COPD

# ------------------------------------------------------------------------------
# "Strong SNPs" score
# ------------------------------------------------------------------------------

#Calculate GRS profile for CAS based on best SNP list

plink --bfile $BFILE --score "$OUTDIR"/sumstats/asthma_grs_weights_strong.txt 5 2 3 sum \
	--memory 64000 --out "$OUTDIR"/asthma.strong

plink --bfile $BFILE --score "$OUTDIR"/sumstats/allergic_grs_weights_strong.txt 4 2 3 sum \
	--memory 64000 --out "$OUTDIR"/allergic.strong
