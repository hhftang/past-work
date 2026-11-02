#!/bin/bash

# ==============================================================================
# GET GRS SCORES FOR SHRINE ET AL SUMMARY STATISTICS    
# ==============================================================================

# rawsummary statistics in 
# /projects/GenomicPrediction/GRS_resources/raw_SummaryStatistics

INDIR=/projects/GenomicPrediction/GRS_resources/COPD_spirometa_ukb
OUTDIR=/home/htang/projects/htang/CAS/C2_genetics/GRS
BFILE=/projects/qinqinhuang/CAS/Genotype_Data/process_data/Michigan_imputation_server_results/cas_auto_sub_filt-updated_auto_chr_2018_01_11_01_15/filtered_geno/cas_michigan_imp_filtered

# /projects/Microbiome/CAS_GRS/ for completed CAS scores from Shumei for selected traits

# ------------------------------------------------------------------------------
# Spirometa traits
# ------------------------------------------------------------------------------

# Format summary statistics file 
# Original format;
# SNP (RSID), Chromosome, Position_b37, Coded (alternate allele), Non-coded (Reference allele), N, Neff, Coded_freq, beta, Se, P
# Desired format
# RSID, alternate allele, beta, Chr:position, Reference allele
gzip -cd "$INDIR"/SpiroMeta_FEV1_results.txt.gz |\
	awk '{ printf $1"\t"$4"\t"$9"\t"$2":"$3"\t"$5"\n"}' > \
	"$OUTDIR"/sumstats/SpiroMeta_FEV1.tsv
# 19814084 rows (including header)

# Take the subset of SNPs that matches UKB r2<0.9 LD-thinned set
grep -w -F -f \
	"$OUTDIR"/ukb_imp_allchr_v2_maf0.001_info0.01_whitebritish_tuning_r2_0.9_noambig.prune.in \
	"$OUTDIR"/sumstats/SpiroMeta_FEV1.tsv > \
	"$OUTDIR"/sumstats/SpiroMeta_FEV1_grs_weights.txt
# 5172191 rows

# Check for duplicate rsID
cat "$OUTDIR"/sumstats/SpiroMeta_FEV1_grs_weights.txt | awk '{print $1}'| sort | uniq -d > \
	"$OUTDIR"/sumstats/SpiroMeta_FEV1_duplicate_SNPS.txt
# Empty - no duplicates

# Calculate GRS profile for CAS based on best SNP list
# Arguments for score represent: file, column of variant ID, allele ID, beta/score 
plink --bfile $BFILE --score "$OUTDIR"/sumstats/SpiroMeta_FEV1_grs_weights.txt 4 2 3 sum \
 --memory 64000 --out "$OUTDIR"/SpiroMeta_FEV1


# ------------------------------------------------------------------------------
# UKB/meta traits
# ------------------------------------------------------------------------------

# Original format;
# SNP (RSID), Chromosome, Position_b37, Coded (alternate allele), Non-coded (Reference allele),INFO, Coded_freq, beta, SE_GC, P_GC
# Desired format
# RSID, alternate allele, beta, Chr:position, Reference allele
gzip -cd "$INDIR"/UKB_FEV1_results.txt.gz |\
	awk '{ printf $1"\t"$4"\t"$8"\t"$2":"$3"\t"$5"\n"}' > \
	"$OUTDIR"/sumstats/UKB_FEV1.tsv
# 19814084 rows (including header)

# Take the subset of SNPs that matches UKB r2<0.9 LD-thinned set
grep -w -F -f \
	"$OUTDIR"/ukb_imp_allchr_v2_maf0.001_info0.01_whitebritish_tuning_r2_0.9_noambig.prune.in \
	"$OUTDIR"/sumstats/UKB_FEV1.tsv > \
	"$OUTDIR"/sumstats/UKB_FEV1_grs_weights.txt
# 5172191 rows

# Check for duplicate rsID
cat "$OUTDIR"/sumstats/UKB_FEV1_grs_weights.txt | awk '{print $1}'| sort | uniq -d > \
	"$OUTDIR"/sumstats/UKB_FEV1_duplicate_SNPS.txt
# Empty - no duplicates

# Calculate GRS profile for CAS based on best SNP list
# Arguments for score represent: file, column of variant ID, allele ID, beta/score 
plink --bfile $BFILE --score "$OUTDIR"/sumstats/UKB_FEV1_grs_weights.txt 4 2 3 sum \
 --memory 64000 --out "$OUTDIR"/UKB_FEV1