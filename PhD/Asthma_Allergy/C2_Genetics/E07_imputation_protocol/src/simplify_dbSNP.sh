#!/bin/bash
# ----------------------
# Simplify snp142.txt from UCSC into something we can work with
# ----------------------
# Keep only those lines that have biallelic, single-base-pair SNPs
# ----------------------
cat snp142.txt | \
  awk -F '\t' '($12 ~ /single/) {print;}' | \
  awk -F '\t' '(length($8)==1) {print;}' > snp142_single_base.txt

# Trim down to columns of interest
# ----------------------
cat snp142_single_base.txt | awk -F '\t' '{print $2, $3, $4, $5, $10}' \
  > dbSNP_database.txt

# Split into chromosomes
# ----------------------
for i in {1..22}; do
	echo "Creating chromosome "$i" file..."
	cat dbSNP_database.txt | egrep "^chr$i[[:space:]]+" | \
	sed 's/chr//' > dbSNP_chr"$i".txt
done

