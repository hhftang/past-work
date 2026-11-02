#!/bin/bash
# ----------------------
# TRIM DOWN TO AUTOSOMES ONLY
# ----------------------
# V. 2015-12-18
# SCRIPT MUST BE IN SAME DIRECTORY AS PLINK FILES.
# Check permissions before running; chmod 744 if necessary.
#
# THIS WORKFLOW TRIMS INPUT PLINK-FORMAT FILES TO INCLUDE ONLY AUTOSOMES
#
# ----------------------
# Set response and instructions if nothing entered
# ----------------------
if [ "$#" -lt 1 ];
then
   echo "auto_only.sh: missing input file"
   echo "usage: ./auto_only.sh <INPUT>"
   exit 1
fi

# Set inputs
# ----------------------
INPUT=$1

# Set other variables
# ----------------------
PLINK="/software/plink-ng/2015-12-13/plink"

# Split file into chromosomes
# Ignore chromosomes "0"(="trashcan"?), "23" (=X), "24" (=Y), "25" (=XY) and 
# "26" (=MT)
# ----------------------
for i in {1..22}; do
	$PLINK --bfile "$INPUT" --chr "$i" --make-bed --out "$INPUT"_chr"$i"
done

# Merge autosomes back into single file
# ----------------------
echo > merge.txt
for i in {2..22}; do
	echo ""$INPUT"_chr"$i".bed "$INPUT"_chr"$i".bim "$INPUT"_chr"$i".fam" \
	  >> merge.txt
done
$PLINK --bfile "$INPUT"_chr1 --merge-list merge.txt --make-bed \
  --out "$INPUT"_auto

# Tidy
# ----------------------
mkdir "$INPUT"_chr
for i in {1..22}; do
	mv "$INPUT"_chr"$i".bed "$INPUT"_chr/.
	mv "$INPUT"_chr"$i".bim "$INPUT"_chr/.
	mv "$INPUT"_chr"$i".fam "$INPUT"_chr/.	
	mv "$INPUT"_chr"$i".log "$INPUT"_chr/.		
done
