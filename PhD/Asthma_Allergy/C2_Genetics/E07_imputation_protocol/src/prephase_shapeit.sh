#!/bin/bash
# ----------------------
# PREIMPUTATION PREPHASING WITH SHAPEIT2 - WORKFLOW
# ----------------------
# V. 2015-12-20
# SCRIPT MUST BE IN SAME DIRECTORY AS PLINK FILES.
# Check permissions before running; chmod 744 if necessary.
#
# THIS WORKFLOW DOES PREPHASING USING SHAPEIT v2.
# CHECKS FOR STRAND ALIGNMENT AND REMOVES DUPLICATE SNPS.
# Make sure data is in hg19/build 37. THIS SCRIPT DOES NOT PERFORM LIFTOVER.
# THIS SCRIPT IS APPROPRIATE FOR AUTOSOMES ONLY
#
# With thanks to Oneil Bhalala and Lesley Raven
# ----------------------
# Set response and instructions if nothing entered
# ----------------------
if [ "$#" -lt 2 ]; then
	echo "prephase_shapeit.sh: missing input file, chromosome"
	echo "usage: ./prephase_shapeit.sh <INPUT> <CHR>"
	echo "Exiting..."
	exit 1
fi

# Set inputs
# ----------------------
INPUT=$1
CHR=$2

# Set other variables
# ----------------------
PLINK="/software/plink-ng/2015-12-13/plink"
SHAPEIT="/software/shapeit.v2.r644.linux.x86_64/shapeit.v2.r644.linux.x86_64"

# Reference data = 1000 Genomes Phase 3 build 37 release Oct 2014
# Check presence in Lesley's folder; exit if absent
# ----------------------
REF="/storage/53145-stor-01/lesley/impute/1000GP_Phase3"

if [ ! -d "$REF" ]; then
	echo "Reference directory missing"
	echo "Consider checking and altering this script as necessary"
	echo "Exiting..."
	exit 1
fi

# Label and extract chromosome in binary plink format
# ----------------------
TAG="$INPUT"_chr"$CHR"
$PLINK --bfile "$INPUT" --chr "$CHR" --make-bed --out "$TAG"

# Check strand alignments
# ----------------------
$SHAPEIT \
  -check \
  -B "$TAG" \
  -M "$REF"/genetic_map_chr"$CHR"_combined_b37.txt \
  -R "$REF"/1000GP_Phase3_chr"$CHR".hap "$REF"/1000GP_Phase3_chr"$CHR".legend \
  "$REF"/1000GP_Phase3.sample \
  -L "$TAG"_alignments \
  -T 30

# If there is an alignment problem:
# Sift out missing and strand alignment errors
# Exclude missing and duplicated snps
# Flip alignment error snps
# ----------------------  
if [ -s "$TAG"_alignments.snp.strand ]; then
	echo "There are alignment issues in Chr "$CHR"; making clean files..."
	awk '{if($1=="missing") print $3}' "$TAG"_alignments.snp.strand | sort > \
	  "$TAG"_missing_snps
	awk '{print $2}' "$TAG".bim | uniq -d | sort > "$TAG"_dup_snps \
	cat "$TAG"_missing_snps "$TAG"_dup_snps > "$TAG"_exclude_snps
	awk '{if($1=="strand") print $3}' "$TAG"_alignments.snp.strand | sort > \
	  "$TAG"_flip_snps
	$PLINK --bfile "$TAG" --exclude "$TAG"_exclude_snps --make-bed \
	  --out "$TAG"_exclude
	$PLINK --bfile "$TAG"_exclude --flip "$TAG"_flip_snps --make-bed --out "$TAG"_clean
else
	echo "No alignment issues in Chr "$CHR"; proceeding as planned..."
	awk '{print $2}' "$TAG".bim | uniq -d | sort > "$TAG"_dup_snps \
	cp "$TAG"_dup_snps "$TAG"_exclude_snps
	$PLINK --bfile "$TAG" --exclude "$TAG"_exclude_snps --make-bed \
	  --out "$TAG"_exclude
	cp "$TAG"_exclude "$TAG"_clean
fi

# Exclude remaining mal-alignments from clean files
# ----------------------
$SHAPEIT \
  -check \
  -B "$TAG"_clean \
  -M "$REF"/genetic_map_chr"$CHR"_combined_b37.txt \
  -R "$REF"/1000GP_Phase3_chr"$CHR".hap "$REF"/1000GP_Phase3_chr"$CHR".legend \
  "$REF"/1000GP_Phase3.sample \
  -L "$TAG"_clean_alignments \
  -T 30

# Prephase *really* clean files
# ----------------------
$SHAPEIT \
  -B "$TAG"_clean \
  -M "$REF"/genetic_map_chr"$CHR"_combined_b37.txt \
  -R "$REF"/1000GP_Phase3_chr"$CHR".hap "$REF"/1000GP_Phase3_chr"$CHR".legend \
  "$REF"/1000GP_Phase3.sample \
  --exclude-snp "$TAG"_clean_alignments.snp.strand.exclude \
  -O "$TAG"_prephased \
  -L "$TAG"_prephased \
  -T 30

# Clean-up
# ----------------------
rm "$TAG"_missing_snps
rm "$TAG"_dup_snps
rm "$TAG"_flip_snps
rm "$TAG"_exclude_snps
# rm "$TAG"_exclude.*
# rm "$TAG"_alignments*
# rm "$TAG"_clean_alignments* 
# rm "$TAG"_clean*

# Tidy-up
# ----------------------
mkdir -p "$INPUT"_shapeit
mv "$TAG".* "$INPUT"_shapeit/
mv "$TAG"_prephased* "$INPUT"_shapeit/