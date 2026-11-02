#!/bin/bash
# ----------------------
# COMBINE CHROMOSOMES
# ----------------------
# Set inputs
# ----------------------
INPUT=$1

# Set other variables
# ----------------------
PLINK="/software/plink-ng/2015-12-13/plink"
OUTDIR="$INPUT"_impute

rm -f mergelist.txt
for i in {1..22}; do
  echo "$OUTDIR"/Chr"$i".combined.{bed,bim,fam} | sed 's/fam /fam\n/g' \
    >> mergelist.txt
done

$PLINK --bfile "$OUTDIR"/Chr22.combined \
   --merge-list mergelist.txt --make-bed --out combined

if [ -s combined-merge.missnp ]
 then
   	echo "There exist missing snps during merge."
    echo "Removing missing snps and then merging..."
   	for((chr=1;chr<=22;chr++))
    	 do
    	 TAG=Chr"$chr"
    	 $PLINK --bfile "$OUTDIR"/"$TAG".combined \
         --exclude combined-merge.missnp \
         --make-bed --out "$OUTDIR"/"$TAG".combined.excl
   	done

   	# echo \
    # 	 "$OUTDIR"/Chr{1,2,3,4,5,6,7,8,9,10,11,\ 
    #      12,13,14,15,16,17,18,19,20,21,22}.combined.excl.{bed,bim,fam} \
    # 	 | sed 's/fam /fam\n/g' > mergelist_excl.txt

    for i in {1..22}; do
      echo "$OUTDIR"/Chr"$i".combined.excl.{bed,bim,fam} |\
        sed 's/fam /fam\n/g' >> mergelist_excl.txt
    done

   	$PLINK --bfile "$OUTDIR"/Chr22.combined.excl \
    	 --merge-list mergelist_excl.txt --make-bed --out combined
   	echo "New combined file made without missing snps."

 else
   	echo "No missing snps present during merge."
fi