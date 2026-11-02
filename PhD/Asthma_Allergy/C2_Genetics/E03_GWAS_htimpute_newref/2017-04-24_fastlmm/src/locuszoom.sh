#!/bin/bash
# ----------------------
# LOCUS ZOOM SCRIPT

TARGET=$1

cat "$TARGET" | cut -f 1,2,4,6 | awk '{print $2,$3,$4,$1}' |\
	cut -f1 -d":" | sort -k 2,2 -g  > locuszoom_output.txt
echo 'Chr BP P.value MarkerName'| cat - locuszoom_output.txt \
	> temp && mv temp locuszoom_output.txt

# Go to http://locuszoom.sph.umich.edu/locuszoom/genform.php?type=yourdata and
# upload locuszoom_output.txt.
#
# Leave all settings default except change to column delimiter to "space"
# Set target SNP +/- 500 kb
# Hg19/1000 Genomes Nov 2014 EUR unless otherwise