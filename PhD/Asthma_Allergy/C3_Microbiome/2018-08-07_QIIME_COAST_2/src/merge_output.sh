#!/bin/bash

#Partition: 'Short' (<24h) or 'Long' (>24h)
#SBATCH -p Short

#Give the job a name:
#SBATCH --job-name="merge"

# Maximum number of tasks/CPU cores used by the job:
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5

# The amount of memory in megabytes per process in the job:
#SBATCH --mem=100000

# The maximum running time of the job in days-hours:mins:sec
#SBATCH --time=0-23:00:00

# Error and output messages
#SBATCH --error=/projects/htang/COAST/2017-12-08_QIIME_COAST/log/merge-output%J.err
#SBATCH --output=/projects/htang/COAST/2017-12-08_QIIME_COAST/log/merge-output%J.out

# ----------------------
# MERGE DADA+FILTERED OUTPUT
# ----------------------

# Merges feature tables & sequences from multiple samples iteratively.

# Set environment
# ----------------------
source activate /usr/local/miniconda3/envs/qiime2-2017.10

# Set temp directory
# ----------------------
# export TMPDIR=/projects/htang/tmp

# Check inputs
# ----------------------

if [ "$#" -lt 1 ];
then
   echo "Missing run list file"
   echo "usage: (sbatch) script.sh <runlist>"
   exit 1
fi

# Set inputs
# ----------------------
src="/projects/htang/COAST/2017-12-08_QIIME_COAST/src/"
data="/projects/htang/COAST/2017-12-08_QIIME_COAST/data/"
runlist=$1

# Combine feature tables
# ----------------------
# Combines a pair of feature tables which contain different samples, and
# which may or may not contain the same features.

# Input: table1 FeatureTable[Frequency]
# Input: table2 FeatureTable[Frequency]
# Output:  merged-table FeatureTable[Frequency]

# Define first table1 as first run

cp $data/AGRF_CAGRF13540_002-filt-table.qza $data/table1.qza

while IFS='' read -r line || [[ -n "$line" ]]; do

  qiime feature-table merge \
    --i-table1 $data/table1.qza \
    --i-table2 $data/$line-filt-table.qza \
    --p-overlap-method sum \
    --o-merged-table $data/table3.qza

  cp $data/table1.qza $data/table1-$line.qza
  mv $data/table3.qza $data/table1.qza

done < $runlist

mv $data/table1.qza $data/table-coast.qza

# Combine feature sequences
# ----------------------
# Combines a pair of feature data objects which may or may not contain data
# for the same features. 

# NOTE: If different feature data is present for the same
# feature id in the two inputs, the data from the first (data1) will be
# propagated to the result.

# Input: data1 FeatureData[Sequence]
# Input: data2 FeatureData[Sequence]
# Output: merged-data FeatureData[Sequence]

cp $data/AGRF_CAGRF13540_002-rep-seqs.qza $data/data1.qza

while IFS='' read -r line || [[ -n "$line" ]]; do

  qiime feature-table merge-seq-data \
    --i-data1 $data/data1.qza \
    --i-data2 $data/$line-rep-seqs.qza \
    --o-merged-data $data/data3.qza

  cp $data/data1.qza $data/data1-$line.qza
  mv $data/data3.qza $data/data1.qza

done < $runlist

mv $data/data1.qza $data/rep-seqs-coast.qza

# Reset temp directory
# ----------------------
export TMPDIR=

# Exit environment
# ----------------------
source deactivate
