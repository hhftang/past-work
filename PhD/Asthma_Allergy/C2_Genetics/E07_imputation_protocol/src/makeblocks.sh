#!/bin/bash
# ----------------------
# MAKE BLOCKS FOR IMPUTATION
# ----------------------
# V. 2016-01-14
# SCRIPT MUST BE IN SAME DIRECTORY AS PLINK FILES.
# Check permissions before running; chmod 744 if necessary.
#
# THIS WORKFLOW MAKES 5MB BLOCKS (NEEDED WITH SHAPEIT PREPHASING) 
# AND CHR-SPECIFIC IMPUTE FILES IN PREPARATION FOR LATER RECONVERSION BACK TO 
# PLINK FORMAT (.SAMPLE FILE FORMAT NEEDED BY GTOOLS)
# 
# With thanks to Oneil Bhalala and Lesley Raven
# ----------------------
# Set response and instructions if nothing entered
# ----------------------
if [ "$#" -lt 1 ]; then
	echo "makeblocks.sh: missing input file"
	echo "usage: ./prephase_shapeit.sh <INPUT> [STEPSIZE]"
	echo "(default stepsize is 5MB)"
	echo "Exiting..."
	exit 1
fi

# Set inputs and stepsize
# ----------------------
INPUT=$1
STEP=$2

if [ -z "$STEP" ]; then
	STEP=5000000
fi

# Set other variables and tools
# ----------------------
GTOOL="/software/gtool_v0.7.5_x86_64/gtool"
PLINK="/software/plink-ng/2015-12-13/plink"

# Convert plink beds to impute format, for each chromosome
# Put all commands into function "run" so that it can be semi-automated
# ----------------------
function run {

# Set chromosome and file tag
# ----------------------
	CHR=$1
	TAG="$INPUT"_"$CHR"

# Identify all lines (SNPs) in input plink file (bim) that belong to the
# specific chromosome. Search for pattern "any_text $CHR space more_text".
# Egrep interprets pattern as extended regular expression.
# Get 4th column (= base position) and keep head and tail. 
# We can then identify base position of first ($START) and last ($END) SNP. 
# ----------------------
	REG="^"$CHR"[[:space:]]+"
	START=$(egrep "$REG" "$INPUT".bim | awk '{print $4}' | head -1)
	END=$(egrep "$REG" "$INPUT".bim | awk '{print $4}' | tail -1)

# Make chromosome-specific PEDs from BEDs
# "2>&1 > /dev/null" means to take all log/error messages and discard them
# ----------------------
	$PLINK --noweb --bfile "$INPUT" --recode --chr "$CHR" --out "$TAG" 2>&1 \
	  > /dev/null
   
# Make impute format from PEDs and get .sample files
# ----------------------
	$GTOOL -P --ped "$TAG".ped --map "$TAG".map \
      --og "$TAG".gen --os "$TAG".sample 2>&1 > /dev/null

# Make "coordinate" file (blocks.txt)
# ----------------------
	while [ "$START" -le "$END" ]; do
      echo "$CHR" "$START" $(($START+$STEP))
      START=$(($START+$STEP))
  	done
}

# Perform function "run" for all chromosomes 1 to 22
# ----------------------
export -f run
export INPUT STEP GTOOL PLINK
seq 1 22 | xargs -I{} bash -c "run {}" > blocks.txt

# Tidy
# ----------------------
mkdir -p "$INPUT"_chr
for i in {1..22}; do
	mv "$INPUT"_"$i".gen "$INPUT"_chr/.
	mv "$INPUT"_"$i".log "$INPUT"_chr/.
	mv "$INPUT"_"$i".map "$INPUT"_chr/.
	mv "$INPUT"_"$i".ped "$INPUT"_cgr/.
	mv "$INPUT"_"$i".sample "$INPUT"_chr/.
done