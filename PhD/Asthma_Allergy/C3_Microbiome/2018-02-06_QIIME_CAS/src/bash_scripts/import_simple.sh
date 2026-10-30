#!/bin/bash

#Partition: 'Short' (<24h) or 'Long' (>24h)
#SBATCH -p Short

#Give the job a name:
#SBATCH --job-name="fqimport"

# Maximum number of tasks/CPU cores used by the job:
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5

# The amount of memory in megabytes per process in the job:
#SBATCH --mem=256000

# The maximum running time of the job in days-hours:mins:sec
#SBATCH --time=0-23:00:00

# Error and output messages
#SBATCH --error=/projects/htang/CAS/C3_microbiome/2017-11-29_QIIME_CAS/log/fqimport%J.err
#SBATCH --output=/projects/htang/CAS/C3_microbiome/2017-11-29_QIIME_CAS/log/fqimport%J.out

# ----------------------
# IMPORT FASTQ FILES TO QIIME2 ARTIFACT FORMAT
# ----------------------

# Set inputs
# ----------------------
src="/projects/htang/CAS/C3_microbiome/2017-11-29_QIIME_CAS/src/"
data="/projects/htang/CAS/C3_microbiome/2017-11-29_QIIME_CAS/data/"
# $1="runs_nogz.txt"
line=$1

# Check inputs
# ----------------------

if [ "$#" -lt 1 ];
then
   echo "Missing input directory"
   echo "usage: (sbatch) script.sh <input>"
   exit 1
fi

# Set environment
# ----------------------
source activate /usr/local/miniconda3/envs/qiime2-2017.10

# Set temp directory
# ----------------------
export TMPDIR=/projects/htang/tmp

# Run import
# ----------------------
# Import data to create a new QIIME 2 Artifact. 
# Target fastq files are CASAVA 1.8 PAIRED-END DEMULTIPLEXED FASTQ files

# Input: FastQ file (CASAVA1.8)
# Output: SampleData[PairedEndSequencesWithQuality] SingleLanePerSamplePairedEndFastqDirFmt

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path $data/$line \
  --source-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path $data/$line-demux-paired-end.qza

# Reset temp directory
# ----------------------
export TMPDIR=

# Exit environment
# ----------------------
source deactivate
