#!/bin/bash

#Partition: 'Short' (<24h) or 'Long' (>24h)
#SBATCH -p Short

#Give the job a name:
#SBATCH --job-name="filthuman"

# Maximum number of tasks/CPU cores used by the job:
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5

# The amount of memory in megabytes per process in the job:
#SBATCH --mem=100000

# The maximum running time of the job in days-hours:mins:sec
#SBATCH --time=0-23:00:00

# Array
#SBATCH --array=1-38%6

# Error and output messages
#SBATCH --error=/projects/htang/CAS/C3_microbiome/2017-11-29_QIIME_CAS/log/filthuman%a_%A.err
#SBATCH --output=/projects/htang/CAS/C3_microbiome/2017-11-29_QIIME_CAS/log/filthuman%a_%A.out

# ----------------------
# QUALITY CONTROL (REMOVE HUMAN READS)
# ----------------------

# This QIIME 2 plugin supports methods for assessing and controlling the quality of feature and sequence data.

# This method gets human reference data, finds sequences in query data
# that matches human data, and excludes those sequences to generate a new
# sequence table


# Set inputs
# ----------------------
src="/projects/htang/CAS/C3_microbiome/2017-11-29_QIIME_CAS/src/"
data="/projects/htang/CAS/C3_microbiome/2017-11-29_QIIME_CAS/data/"
ref="/projects/htang/human_ncbi_fasta/human_ref.qza"
runlist=$1
query=`head -n $SLURM_ARRAY_TASK_ID $runlist | tail -1`

# Check inputs
# ----------------------

if [ "$#" -lt 1 ];
then
   echo "Missing runlist"
   echo "usage: (sbatch) script.sh <runlist>"
   exit 1
fi

# Go to directory
# ----------------------
# cd /projects/htang/CAS/C3_microbiome/

# Set environment
# ----------------------
source activate /usr/local/miniconda3/envs/qiime2-2017.10

# Set temp directory
# ----------------------
export TMPDIR=/projects/htang/tmp

# Get human sequences to exclude (hits)
# ----------------------
# This method aligns feature sequences to a set of reference sequences to
# identify sequences that hit/miss the reference within a specified
# perc_identity, evalue, and perc_query_aligned. This method could be used
# to define a positive filter, e.g., extract only feature sequences that
# align to a certain clade of bacteria; or to define a negative filter,
# e.g., identify sequences that align to contaminant or human DNA sequences
# that should be excluded from subsequent analyses. Note that filtering is
# performed based on the perc_identity, perc_query_aligned, and evalue
# thresholds (the latter only if method==BLAST and an evalue is set). Set
# perc_identity==0 and/or perc_query_aligned==0 to disable these filtering
# thresholds as necessary.

# Input: query-sequences: FeatureData[Sequence] 
# Input: reference-sequences: FeatureData[Sequence] 
# Output: sequence-hits: FeatureData[Sequence] 
# Output: sequence-misses: FeatureData[Sequence]  

qiime quality-control exclude-seqs \
  --i-query-sequences $data/$query-rep-seqs.qza \
  --i-reference-sequences $ref \
  --p-method blast \
  --p-perc-identity 0.97 \
  --p-perc-query-aligned 0.97 \
  --o-sequence-hits $data/$query-hits.qza \
  --o-sequence-misses $data/$query-misses.qza

# Get human sequences to exclude (hits)
# ----------------------
# Filter features from table based on frequency and/or metadata. Any samples
# with a frequency of zero after feature filtering will also be removed. See
# the filtering tutorial on https://docs.qiime2.org for additional details.

# Input: table: FeatureTable[Frequency]
# Output: filtered-table: FeatureTable[Frequency]

qiime feature-table filter-features \
  --i-table $data/$query-dada-table.qza  \
  --m-metadata-file $data/$query-hits.qza \
  --o-filtered-table $data/$query-filt-table.qza \
  --p-exclude-ids

# Reset temp directory
# ----------------------
export TMPDIR=

# Exit environment
# ----------------------
source deactivate
