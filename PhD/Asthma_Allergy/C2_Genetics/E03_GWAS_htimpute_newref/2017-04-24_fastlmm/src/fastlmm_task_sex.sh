#!/bin/bash
# ----------------------
# GWAS USING FASTLMM WITH SEX AS COVARIATE
# ----------------------
# V. 2016-08-09
# Modified to include sex
# SCRIPT MUST BE IN SAME DIRECTORY AS PLINK FILES.
#
# PLACE PHENO FILES IN /PHENO/. DIRECTORY.
# ----------------------
# Set response and instructions if nothing entered
# ----------------------
if [ "$#" -lt 2 ]; then
	echo "fastlmm_task_sex.sh: missing input file, task"
	echo "usage: ./fastlmm_task_sex.sh <INPUT> <TASK>"
	echo "Exiting..."
	exit 1
fi

# Set inputs
# ----------------------
INPUT=$1
TASK=$2

# Set other variables
# ----------------------
PLINK="/software/plink-ng/2015-12-13/plink"
FASTLMMC="/software/FaSTLMM.207/Bin/Linux_MKL/fastlmmc"

# Perform fastlmm
# ----------------------
for i in {1..22}; do
	mkdir -p output/"$TASK"_sex
	$FASTLMMC -verboseOutput -maxThreads 20 \
	  -bfile "$INPUT"_Chr"$i" -bfilesim "$INPUT"_nChr"$i" -pheno pheno/"$TASK".phe \
	  -covar pheno/sex.phe -out output/"$TASK"_sex/"$TASK"_sex_Chr"$i".txt

# Trim out headers from final file
# ----------------------
	mv output/"$TASK"_sex/"$TASK"_sex_Chr"$i".txt output/"$TASK"_sex/"$TASK"_sex_old_Chr"$i".txt
	tail -n+2 output/"$TASK"_sex/"$TASK"_sex_old_Chr"$i".txt > output/"$TASK"_sex/"$TASK"_sex_Chr"$i".txt
done

# Grab header of original file format, and clean-up
# ----------------------
head -1 output/"$TASK"_sex/"$TASK"_sex_old_Chr1.txt > output/"$TASK"_sex/"$TASK"_sex_old_header.txt
rm output/"$TASK"_sex/"$TASK"_sex_old_Chr*.txt 

# Combine files to form master file
# ----------------------
cat output/"$TASK"_sex/"$TASK"_sex_Chr*.txt | cut -f 1,2,4,6 > output/"$TASK"_sex/master_"$TASK"_sex.txt
sort -gk4,4 output/"$TASK"_sex/master_"$TASK"_sex.txt > output/"$TASK"_sex/master_"$TASK"_sex_sort.txt

