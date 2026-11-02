#!/bin/bash

# ==============================================================================
# GET GRS SCORES FOR SHRINE ET AL SUMMARY STATISTICS - Spirometa
# ==============================================================================

# rawsummary statistics in 
# /projects/GenomicPrediction/GRS_resources/raw_SummaryStatistics

PHENO=$1
INDIR=/projects/GenomicPrediction/GRS_resources/COPD_spirometa_ukb
OUTDIR=/projects/htang/CAS/C2_genetics/GRS
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
gzip -cd "$INDIR"/"$PHENO"_results.txt.gz |\
	awk '{ printf $1"\t"$4"\t"$9"\t"$2":"$3"\t"$5"\n"}' > \
	"$OUTDIR"/sumstats/"$PHENO".tsv

# Take the subset of SNPs that matches UKB r2<0.9 LD-thinned set
grep -w -F -f \
	"$OUTDIR"/ukb_imp_allchr_v2_maf0.001_info0.01_whitebritish_tuning_r2_0.9_noambig.prune.in \
	"$OUTDIR"/sumstats/"$PHENO".tsv > \
	"$OUTDIR"/sumstats/"$PHENO"_grs_weights.txt

# Check for duplicate rsID
cat "$OUTDIR"/sumstats/"$PHENO"_grs_weights.txt | awk '{print $1}'| sort | uniq -d > \
	"$OUTDIR"/sumstats/"$PHENO"_duplicate_SNPS.txt

# Calculate GRS profile for CAS based on best SNP list
# Arguments for score represent: file, column of variant ID, allele ID, beta/score 
plink --bfile $BFILE --score "$OUTDIR"/sumstats/"$PHENO"_grs_weights.txt 4 2 3 sum \
 --memory 64000 --out "$OUTDIR"/"$PHENO"

