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
#SBATCH --error=/projects/htang/CAS/C3_microbiome/2017-11-29_QIIME_CAS/log/vsearch%J.err
#SBATCH --output=/projects/htang/CAS/C3_microbiome/2017-11-29_QIIME_CAS/log/vsearch%J.out

# ----------------------
# Get closed reference OTU picking so that we can actually do copy number correction...
# ----------------------

# Set environment
# ----------------------
source activate /usr/local/miniconda3/envs/qiime2-2017.10

# Set temp directory
# ----------------------
export TMPDIR=/projects/htang/tmp

# Set inputs
# ----------------------
data="/projects/htang/CAS/C3_microbiome/2017-11-29_QIIME_CAS/data/"

# Steps done before hand
# qiime tools import \
# 	--input-path "$data"/gg_13_5_otus/rep_set/99_otus.fasta \
# 	--output-path "$data"/gg_13_5_otu_99.qza \
# 	--type 'FeatureData[Sequence]'

# qiime feature-table tabulate-seqs \
#   --i-data "$data"/gg_13_5_otu_99.qza \
#   --o-visualization "$data"/gg_13_5_otu_99.qzv

# NOTE: THIS PARICULAR STEP REQUIRES QIIME2-2017.12
# qiime feature-table filter-seqs \
# 	--i-data "$data"/rep-seqs-cas.qza \
# 	--i-table "$data"/table-cas.qza \
# 	--o-filtered-data "$data"/rep-seqs-cas-filt.qza

qiime vsearch cluster-features-closed-reference \
	--i-sequences "$data"/rep-seqs-cas-filt.qza \
	--i-table "$data"/table-cas.qza \
	--i-reference-sequences "$data"/gg_13_5_otu_99.qza \
	--p-perc-identity 1 \
	--p-threads 0 \
	--output-dir "$data"/closedref

# Reset temp directory
# ----------------------
export TMPDIR=

# Exit environment
# ----------------------
source deactivate
