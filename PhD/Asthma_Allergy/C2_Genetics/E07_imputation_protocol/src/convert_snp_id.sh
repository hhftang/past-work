#!/bin/bash
# ----------------------
# CONVERT KGP IDS TO RS IDS
# ----------------------
# V. 2015-12-20
# SCRIPT MUST BE IN SAME DIRECTORY AS PLINK FILES.
# Check permissions before running; chmod 744 if necessary.
#
# THIS WORKFLOW DOES THE FOLLOWING:
# 1. REMOVES KGP (ILLUMINA) SNPS THAT ARE DUPLICATES OF RS (DBSNP) SNPS
# 2. CONVERTS AS MANY OF THE REMAINING KGP SNPS TO RS SNP IDS 
#
# NOTE: MAKE SURE YOU INCLUDE CONVERSION FILE IN SAME DIRECTORY AS PLINK FILES
# e.g. illumina_kgp_to_rsID.txt
# ----------------------
# Set response and instructions if nothing entered
# ----------------------
if [ "$#" -lt 1 ];
then
   echo "convert_snp_id.sh: missing input file"
   echo "usage: convert_snp_id.sh <INPUT>"
   exit 1
fi

# Set inputs
# ----------------------
INPUT=$1

# Set other variables
# ----------------------
PLINK="/software/plink-ng/2015-12-13/plink"
DBSNP="/storage/abyss/hhtang-u/UCSC_hg19_dbSNP"

# 1. REMOVE DUPLICATES
# Identify SNP duplicates (same chromosome number, same position) with kgp or 
# alternate Illumina IDs.
#
# For 1st bim file (denoted by NR==FNR i.e. total number of records that have 
# been read must match current record number; FNR resets once new file is read): 
# Create an array "a" with dimensions Chr ($1) and Chr pos ($4). 
# For each instance of $1,$4, add one to that cell a[$1,$4].
#
# For 2nd bim file:
# If cell a[$1,$4]>1 (i.e. duplicate exists), where snp id ($2) contains kgp or
# VG, print $2 into file for exclusion.
# ----------------------
awk 'NR==FNR{a[$1,$4]++;next;} \
  {if (a[$1,$4] >1 && ($2 ~/kgp/ || $2 ~/VG/)) print $2;}' \
  "$INPUT".bim "$INPUT".bim > kgp_dup_to_exclude.txt

# Exclude these replicates to create trimmed files
# ----------------------
$PLINK --bfile "$INPUT" --exclude kgp_dup_to_exclude.txt --make-bed \
  --out "$INPUT"_trim

# 2. CONVERSION BY PROVIDED MAPPING FILE
# Use conversion file provided by illumina to map kgp IDs to rs 
# Conversion file name: illumina_kgp_to_rsID.txt
# Note that original file has been edited to remove blank (kgpID to .) 
# mappings and rsID re-mappings
# ----------------------
$PLINK --bfile "$INPUT"_trim --update-name illumina_kgp_to_rsID.txt --make-bed \
  --out "$INPUT"_convert

# 3. CREATION OF MAPPING FILE FOR REMAINING UNSUBBED SNPS
# Find remaining unsubstituted kgps, and identify rsID using existing
# dbSNP database
# ----------------------
grep 'kgp\|VG' "$INPUT"_convert.bim | awk '{print $1, $4-1, $4, $2, $6}' \
  > kgp_unsubbed.txt

# For each chromosome, print kgp-to-rsID matches
# Match by chromosome number, base positions, +/- one of the alleles
# ----------------------
for i in {1..22}; do
	cat kgp_unsubbed.txt | awk -v num=$i '{if($1==num) print;}' \
	  > kgp_unsubbed_chr"$i".txt
done

# Code if we insist on matching alleles:
# For 1st file (kgp list), create two arrays "a" and "b", each with dimensions 
# Chr number ($1), and Chr base positions ($2,$3). 
# ----------------------
# echo > kgp_rs_matching.txt
# for i in {1..22}; do
#   echo "Subbing chromosome ""$i"
#   awk 'NR==FNR{a[$1,$2,$3]=$4; b[$1,$2,$3]=$5; next;} \
# 	  {if((($1,$2,$3) in a) && ($5 ~ b[$1,$2,$3])) 
#     print a[$1,$2,$3], $4, b[$1,$2,$3], $5;}' \
# 	  kgp_unsubbed_chr"$i".txt "$DBSNP"/dbSNP_chr"$i".txt \
# 	  >> kgp_rs_matching.txt
# done

# Code if we do NOT insist on matching alleles:
# ----------------------
echo > kgp_rs_matching.txt
for i in {1..22}; do
  echo "Subbing chromosome ""$i"  
  awk 'NR==FNR{a[$1,$2,$3]=$4; b[$1,$2,$3]=$5; next;} \
    {if(($1,$2,$3) in a) 
    print a[$1,$2,$3], $4;}' \
    kgp_unsubbed_chr"$i".txt "$DBSNP"/dbSNP_chr"$i".txt \
    >> kgp_rs_matching.txt
done
#   print a[$1,$2,$3], $4, b[$1,$2,$3], $5;}' \ # to print alleles also

# Remove any duplicate assignments 
# i.e. if a kgp ID has already been mapped to a rs ID, remove any subsequent 
# assignments with other rsIDs 
# ----------------------
awk '!a[$1]++' kgp_rs_matching.txt > kgp_rs_matching2.txt
mv kgp_rs_matching2.txt kgp_rs_matching.txt

# Convert remaining SNPs using new conversion file:
# ----------------------
$PLINK --bfile "$INPUT"_convert --update-name kgp_rs_matching.txt --make-bed \
  --out "$INPUT"_convert2

# Tidy kgp matching files
# ----------------------
mkdir kgp_rs
mv kgp_unsubbed* kgp_rs/.
mv kgp_rs_matching.txt kgp_rs/.
mkdir "$INPUT"_trim
mv "$INPUT"_trim.* "$INPUT"_trim/.
mv kgp_dup_to_exclude.txt "$INPUT"_trim/.
mkdir "$INPUT"_convert
mv "$INPUT"_convert.* "$INPUT"_convert/.
