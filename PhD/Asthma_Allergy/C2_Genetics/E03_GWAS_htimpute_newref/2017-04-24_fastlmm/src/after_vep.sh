#!/bin/bash

TASK=$1

cut -f1,4,5,6,8  snp_"$TASK"_results.txt > snp_"$TASK"_results2.txt
join -1 1 -2 1 -o 0,1.2,1.3,1.4,1.5,2.2,2.3,2.4 -e NA -t $'\t' \
  <(sort -k1,1 snp_"$TASK"_results2.txt) <(sort -k1,1 master_"$TASK"_fixed.txt) \
 | sort -gk8,8 > snp_"$TASK"_results3.txt

