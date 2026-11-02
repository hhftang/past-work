#!/bin/bash

# Set response and instructions if nothing entered
# ----------------------
if [ "$#" -lt 2 ]; then
	echo "vep.sh: missing task list"
	echo "usage: ./vep.sh <pheno_list>"
	echo "Exiting..."
	exit 1
fi

# Set inputs
# ----------------------
pheno_list=( "$@" )

for i in "${pheno_list[@]}"
do
	paste <(cut -f1 output/"$i"/master_"$i"_sort.txt | cut -f1 -d":") \
  <(cut -f 2,3,4,5 output/"$i"/master_"$i"_sort.txt) \
  | awk '$4<0.00001{print}' > output/"$i"/master_"$i"_fixed.txt
  cut -f1  output/"$i"/master_"$i"_fixed.txt > output/"$i"/snp_"$i".txt

  vep --cahce -i output/"$i"/snp_"$i".txt -o

 	paste <(cut -f1 output/"$i"_sex/master_"$i"_sex_sort.txt | cut -f1 -d":") \
  <(cut -f 2,3,4,5 output/"$i"_sex/master_"$i"_sex_sort.txt) \
  | awk '$4<0.00001{print}' > output/"$i"_sex/master_"$i"_sex_fixed.txt
  cut -f1  output/"$i"_sex/master_"$i"_sex_fixed.txt > output/"$i"_sex/snp_"$i"_sex.txt
 done

 