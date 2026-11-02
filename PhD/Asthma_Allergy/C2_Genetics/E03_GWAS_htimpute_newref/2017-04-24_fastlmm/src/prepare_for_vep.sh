#!/bin/bash

TASK=$1

paste <(cut -f1 master_"$TASK"_sort.txt | cut -f1 -d":") \
  <(cut -f 2,3,4,5 master_"$TASK"_sort.txt) \
  | awk '$4<0.00001{print}' > master_"$TASK"_fixed.txt
cut -f1  master_"$TASK"_fixed.txt > snp_"$TASK".txt