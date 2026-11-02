#!/bin/bash

#Partition: 'Short' (<24h) or 'Long' (>24h)
#SBATCH -p Short

#Give the job a name:
#SBATCH --job-name="sum.out"

# Maximum number of tasks/CPU cores used by the job:
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5

# The amount of memory in megabytes per process in the job:
#SBATCH --mem=256000

# The maximum running time of the job in days-hours:mins:sec
#SBATCH --time=0-23:00:00

# Error and output messages
#SBATCH --error=/projects/htang/COAST/2017-12-08_QIIME_COAST/log/summarize-output%J.err
#SBATCH --output=/projects/htang/COAST/2017-12-08_QIIME_COAST/log/summarize-output%J.out

# ----------------------
# GENERATE SUMMARY VISUALISATION OF QIIME2 SEQUENCE ARTIFACT
# ----------------------

# Set environment
# ----------------------
source activate /usr/local/miniconda3/envs/qiime2-2017.10

# Set temp directory
# ----------------------
export TMPDIR=/projects/htang/tmp

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

# Create combined metadata file
# ----------------------

cp $data/AGRF_CAGRF13540_002.tsv $data/metadata.tsv

while IFS='' read -r line || [[ -n "$line" ]]; do

  tail -n +2 $data/$line.tsv > $data/metadata2.tsv 
	cat $data/metadata.tsv $data/metadata2.tsv > $data/metadata3.tsv
	cp $data/metadata3.tsv $data/metadata.tsv

done < $runlist

rm $data/metadata2.tsv

# Create visualisations
# ----------------------

qiime feature-table summarize \
  --i-table $data/table-coast.qza \
  --o-visualization $data/table-coast.qzv \
  --m-sample-metadata-file $data/metadata.tsv

qiime feature-table tabulate-seqs \
  --i-data $data/rep-seqs-coast.qza \
  --o-visualization $data/rep-seqs-coast.qzv

# Reset temp directory
# ----------------------
export TMPDIR=

# Exit environment
# ----------------------
source deactivate
