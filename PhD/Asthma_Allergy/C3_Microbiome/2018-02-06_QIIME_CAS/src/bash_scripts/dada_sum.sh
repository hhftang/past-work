#!/bin/bash

#Partition: 'Short' (<24h) or 'Long' (>24h)
#SBATCH -p Short

#Give the job a name:
#SBATCH --job-name="cas.dasum"

# Maximum number of tasks/CPU cores used by the job:
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5

# The amount of memory in megabytes per process in the job:
#SBATCH --mem=256000

# The maximum running time of the job in days-hours:mins:sec
#SBATCH --time=0-23:00:00

# ----------------------
# DADA2 SUMMARIZE ONLY
# ----------------------

# This QIIME 2 plugin wraps DADA2 and supports sequence quality control for 
# single-end and paired-end reads using the DADA2 R library.

# This method denoises paired-end sequences, dereplicates them, and filters
# chimeras.


# Set inputs
# ----------------------
INPUT=$1
TRIMLENGTH=$2 # 10
TRUNCLENGTH=$3 # 150

# Check inputs
# ----------------------

if [ "$#" -lt 1 ];
then
   echo "Missing input location, trim length or truncation length"
   echo "usage: (sbatch) script.sh <INPUT> <TRIMLENGTH> <TRUNCLENGTH>"
   exit 1
fi

# Go to directory
# ----------------------
# cd /projects/htang/CAS/C3_microbiome/

# Set environment
# ----------------------
source activate /usr/local/miniconda3/envs/qiime2-2017.9

# Set temp directory
# ----------------------
# export TMPDIR=/projects/htang/tmp

# Run dada2 denoise
# ----------------------

qiime feature-table summarize \
  --i-table "$INPUT"-dada-table.qza \
  --o-visualization "$INPUT"-dada-table.qzv \
  --m-sample-metadata-file "$INPUT"_metadata/SampleSheet_"$INPUT".tsv

# Reset temp directory
# ----------------------
# export TMPDIR=

# Exit environment
# ----------------------
source deactivate
