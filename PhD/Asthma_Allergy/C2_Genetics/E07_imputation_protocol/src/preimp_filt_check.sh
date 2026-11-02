#!/bin/bash
# ----------------------
# PREIMPUTATION FILTER AND CHECK - WORKFLOW
# ----------------------
# V. 2015-12-20
# SCRIPT MUST BE IN SAME DIRECTORY AS PLINK FILES.
# Check permissions before running; chmod 744 if necessary.
#
# THIS WORKFLOW DOES THE FOLLOWING:
# 1. PREIMPUTATION FILTERING (GENO, MIND, MAF, HWE),
# 2. IBD CHECK, HETEROZYGOSITY RATE AND SEX CHECK.
#
# THIS WORKFLOW DOES NOT DO LIFTOVER OR CONVERSION TO FORWARD STRAND.
# IT IS ASSUMED THAT EVERYTHING IS IN THE CORRECT BUILD (37)
# ----------------------
# Set response and instructions if nothing entered
# ----------------------
if [ "$#" -lt 1 ];
then
   echo "preimp_filt_check.sh: missing input file"
   echo "usage: ./preimp_filt_check.sh <INPUT>"
   echo "also remember to have conversion files in /kgprs/chrX_kgprs.txt"
   echo "open and read script for more details"
   exit 1
fi

# Set inputs
# ----------------------
INPUT=$1

# Set other variables
# ----------------------
PLINK="/software/plink-ng/2015-12-13/plink"

# 1. Filter file
# Criteria: include only geno <0.01, mind <0.01, maf >0.01; HWE>1e-6
# ----------------------
$PLINK --bfile "$INPUT" --geno 0.01 --mind 0.01 --maf 0.01 --hwe 0.000001 \
  --make-bed --out "$INPUT"_filt

# Perform LD pruning for IBD analysis
# Generates a subset of SNPs in relative linkage equilibrium with each other
# ----------------------
$PLINK --bfile "$INPUT"_filt --indep-pairwise 1000 50 0.05 --out "$INPUT"_indep
$PLINK --bfile "$INPUT"_filt --extract "$INPUT"_indep.prune.in --make-bed \
  --out "$INPUT"_prune 

# Pairwise IBD estimation, and pi-hat from pruned dataset
# Identify those above pi-hat threshold 0.01
# ----------------------
$PLINK --bfile "$INPUT"_prune --genome --out "$INPUT"_ibd
awk '$10 > 0.125' "$INPUT"_ibd.genome > "$INPUT"_ibd_related.txt

# Hardy-Weinberg and homozygosity statistics
# ----------------------
$PLINK --bfile "$INPUT"_filt --hardy --out "$INPUT"_hwe
$PLINK --bfile "$INPUT"_filt --homozyg --out "$INPUT"_homozyg

# Tidy
# ----------------------
mkdir -p "$INPUT"_check
mv "$INPUT"_homozyg* "$INPUT"_check/
mv "$INPUT"_hwe* "$INPUT"_check/
mv "$INPUT"_ibd* "$INPUT"_check/
mv "$INPUT"_indep* "$INPUT"_check/
mv "$INPUT"_prune* "$INPUT"_check/
