#!/bin/bash
# ----------------------
# POST-IMPUTATION FILTERING AND COMBINE BLOCKS (CHROMOSOME-SPECIFIC)
# ----------------------
# V. 2016-01-15
#
# THIS WORKFLOW WILL CONVERT .IMPUTE FILES, FILTER EACH IMPUTED CHROMOSOME
# AND COMBINE ALL IMPUTED BLOCKS.
# FILTERING REMOVES SNPS WITH MULTIPLE ALLELES, AND LOW POST-IMPUTATION 
# PROBABILITY.
# 
# Imputed files should be in "$INPUT"_impute/.
# Chromosome sample files from makeblocks.sh should be in "$INPUT"_chr/.
#
# With thanks to Oneil Bhalala and Lesley Raven
# ----------------------
# Set inputs
# ----------------------
INPUT=$1
STARTCHR=$2
ENDCHR=$3

# Set other variables
# ----------------------
GTOOL="/storage/abyss/software/gtool_v0.7.5_x86_64/gtool"
PLINK="/software/plink-ng/2015-12-13/plink"
OUTDIR="$INPUT"_impute

for((chr=$STARTCHR;chr<=$ENDCHR;chr++))
do
   TAG=Chr"$chr"
   GEN="$TAG"_combined.gen
   N=$(tail -n +3 "$INPUT"_chr/"$INPUT"_"$chr".sample | wc -l)
   if [ -s "$OUTDIR/$GEN.ped" ] && \
     [ $(wc -l "$OUTDIR/$GEN.ped" | cut -f1 -d' ') -eq $N ]
   then
      echo "skipping $GEN.ped, already exists"
   else
# Concatenate all impute2 and info files from the one chromosome
# impute2 files concatenated into "$TAG"_combined.gen
# info files concatenated into "$TAG"_combined.info
# ----------------------
      cat "$OUTDIR"/"$TAG"_*_*.impute2 > "$OUTDIR"/"$GEN"
      cat "$OUTDIR"/"$TAG"_*_*.impute2_info > "$OUTDIR"/"$TAG"_combined.info

# Include only the mononucleotide alleles (i.e. exclude multiple alleles, and
# polynucleotide SNPs)
# ----------------------
      awk '{ if ($2~"rs." && length($4)==1 && length($5)==1) print }' \
        "$OUTDIR"/"$GEN" > "$OUTDIR"/"$GEN"_2
      awk '{ if ($2~"rs.") print }' "$OUTDIR"/"$TAG"_combined.info > \
        "$OUTDIR"/"$TAG"_combined.info_2

# Quick clean-up
# ----------------------
      rm "$OUTDIR"/"$GEN"
      rm "$OUTDIR"/"$TAG"_combined.info

# Enforce post-imputation filter (info>0.4)
# ----------------------
      awk '{if(NR > 1 && $7 > 0.4) {print $2}}' \
        "$OUTDIR"/"$TAG"_combined.info_2 > \
        "$OUTDIR"/"$TAG"_combined_goodsnps.txt

# Convert impute/gen file format back to ped.
# Probability threshold for assignment of allele is 0.9
# ----------------------
   $GTOOL -G --g "$OUTDIR"/"$GEN"_2 \
     --s "$INPUT"_chr/"$INPUT"_"$chr".sample \
	  --chr "$chr" --snp --threshold 0.9 \
	  --log /dev/null 
   fi

# Put all the .ped together for that chromosome   
# ----------------------
   $PLINK --noweb --file "$OUTDIR"/"$GEN"_2 \
      --missing-genotype N \
      --extract "$OUTDIR"/"$TAG"_combined_goodsnps.txt --make-bed \
      --out "$OUTDIR"/"$TAG".combined

# Clean up
# ----------------------
   rm "$OUTDIR"/"$GEN"_2*
   rm "$OUTDIR"/"$TAG"_combined.info_2
done

#echo \
#   $OUTDIR/chr{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22}.combined.{bed,bim,fam} \
#   | sed 's/fam /fam\n/g' > mergelist_WTCCC.txt

#eval $PLINK --noweb --bfile $OUTDIR/$f.combined \
#   --merge-list mergelist_WTCCC.txt --make-bed --out combined_WTCCC

#eval $PLINK --noweb --bfile combined \
#   --allow-no-sex \
#   --maf 0.01 --hwe 5e-6 --mind 0.05 --geno 0.01 \
#   --make-bed --out combined_filtered

