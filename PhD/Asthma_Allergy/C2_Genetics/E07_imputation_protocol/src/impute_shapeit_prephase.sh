#!/bin/bash
# ----------------------
# IMPUTE BLOCKS WITH PREPHASED FILES
# ----------------------
# V. 2016-01-14
# SCRIPT MUST BE IN SAME DIRECTORY AS BLOCK.TXT FILE.
# Also, make sure shapeit haplotype files are in:
	# $INPUT_shapeit/$INPUT_chr$N_prephased.haps
# and note $INPUT has to be entered manually into the run function below.
# Change SHAPEHAP variable in run function below if necessary.
# Check permissions before running; chmod 744 if necessary.
#
# THIS WORKFLOW IMPUTES THE 5MB BLOCKS CREATED BY MAKEBLOCKS.SH, 
# USING THE SHAPEIT FILES CREATED BY PREPHASE_SHAPEIT.SH
#
# With thanks to Oneil Bhalala and Lesley Raven
# ----------------------
# Immediately stop script if there is an error
# ----------------------
# set -e

# Set response and instructions if nothing entered
# ----------------------
if [ "$#" -lt 1 ]; then
   echo "impute_shapeit_prephase.sh: missing input file or target chromosome"
   echo "usage: ./impute_shapeit_prephase.sh <INPUT> <CHR> [NUMPROCS]"
   exit 1
fi

# Set inputs
# ----------------------
INPUT=$1
CHR=$2
NUMPROCS=$3

# If number of processors is not set, default=1
# ----------------------
if [ -z "$NUMPROCS" ]; then
	NUMPROCS=1
fi

# Set other variables
# ----------------------
IMPUTE="/home/UNIMELB/lraven-u/bin/impute_v2.3.2_x86_64_static/impute2" 

# Reference data = 1000 Genomes Phase 3 build 37 release Oct 2014; in Lesley's folder
# Check presence and exit if absent
# ----------------------
REF="/storage/53145-stor-01/lesley/impute/1000GP_Phase3"

if [ ! -d "$REF" ]; then
	echo "Reference directory missing"
	echo "Consider checking and altering this script as necessary"
	echo "Exiting..."
	exit 1
fi

# Tidy (beforehand)
# ----------------------
OUTDIR="$INPUT"_impute
mkdir -p "$OUTDIR"

# Only identify lines in block.txt with desired chromosome
# For awk, you need to assign shell variable ("$CHR") to a pattern ("pat") 
# that awk can use.
# ----------------------
awk -v pat="$CHR" '{if($1==pat) print}' blocks.txt > blocks_Chr"$CHR".txt

# IMPUTE2 WITH SHAPEIT PREPHASED FILES
# ----------------------
function run {
	BLOCKCHR=$1
	START=$2
	END=$3
	INPUT="cas_auto_sub_filt"	
	HAPLO="$REF"/1000GP_Phase3_chr"$BLOCKCHR".hap
	SHAPEHAP="$INPUT"_shapeit/"$INPUT"_chr"$BLOCKCHR"_prephased.haps
	echo "Chromosome: "$BLOCKCHR"; Block start: "$START""
	echo "Reference haplotype: "$HAPLO""
	echo "Shapeit prephased haplotype file: "$SHAPEHAP""
	
	$IMPUTE \
	  -use_prephased_g \
	  -known_haps_g "$SHAPEHAP"\
	  -m "$REF"/genetic_map_chr"$BLOCKCHR"_combined_b37.txt \
	  -h "$HAPLO" \
	  -l "$REF"/1000GP_Phase3_chr"$BLOCKCHR".legend \
	  -int "$START" "$END" \
	  -align_by_maf_g \
	  -Ne 20000 \
	  -k 100 \
	  -buffer 1000 \
	  -o "$OUTDIR"/Chr"$BLOCKCHR"_"$START"_"$END".impute2
}

# Perform function "run" only for target chromosome, for each block
# ----------------------
export -f run
export CHR NUMPROCS REF OUTDIR IMPUTE
cat blocks_Chr"$CHR".txt | xargs -P"$NUMPROCS" -I{} bash -c "run {}"