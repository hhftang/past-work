#!/bin/bash

#Partition: 'Short' (<24h) or 'Long' (>24h)
#SBATCH -p Short

#Give the job a name:
#SBATCH --job-name="dada"

# Maximum number of tasks/CPU cores used by the job:
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5

# The amount of memory in megabytes per process in the job:
#SBATCH --mem=256000

# The maximum running time of the job in days-hours:mins:sec
#SBATCH --time=0-23:00:00

# Error and output messages
#SBATCH --error=/projects/htang/CAS/C3_microbiome/2017-11-29_QIIME_CAS/log/dada%J.err
#SBATCH --output=/projects/htang/CAS/C3_microbiome/2017-11-29_QIIME_CAS/log/dada%J.out

# ----------------------
# DADA2 DENOISE OF QIIME2 SEQUENCE ARTIFACT
# ----------------------

# This QIIME 2 plugin wraps DADA2 and supports sequence quality control for 
# single-end and paired-end reads using the DADA2 R library.

# This method denoises paired-end sequences, dereplicates them, and filters
# chimeras.


# Set inputs
# ----------------------
src="/projects/htang/CAS/C3_microbiome/2017-11-29_QIIME_CAS/src/"
data="/projects/htang/CAS/C3_microbiome/2017-11-29_QIIME_CAS/data/"
line=$1
trimlength=$2 # 10
trunclength=$3 # 150

# Check inputs
# ----------------------

if [ "$#" -lt 1 ];
then
   echo "Missing input location, trim length or truncation length"
   echo "usage: (sbatch) script.sh <runlist> <trimlength> <trunclength>"
   exit 1
fi

# Set environment
# ----------------------
source activate /usr/local/miniconda3/envs/qiime2-2017.10

# Set temp directory
# ----------------------
export TMPDIR=/projects/htang/tmp

# Run dada2 denoise
# ----------------------
# Input: SampleData[PairedEndSequencesWithQuality]
# Output: dada-table = FeatureTable[Frequency] BIOMV210DirFmt
# Output: rep-seqs = FeatureData[Sequence] DNASequencesDirectoryFormat

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs $data/$line-demux-paired-end.qza \
  --o-table $data/$line-dada-table \
  --o-representative-sequences $data/$line-rep-seqs \
  --p-trim-left-f $trimlength \
  --p-trim-left-r $trimlength \
  --p-trunc-len-f $trunclength \
  --p-trunc-len-r $trunclength

qiime feature-table summarize \
  --i-table $data/$line-dada-table.qza \
  --o-visualization $data/$line-dada-table.qzv \
  --m-sample-metadata-file $data/$line.tsv

qiime feature-table tabulate-seqs \
  --i-data $data/$line-rep-seqs.qza \
  --o-visualization $data/$line-rep-seqs.qzv

# Reset temp directory
# ----------------------
export TMPDIR=

# Exit environment
# ----------------------
source deactivate
