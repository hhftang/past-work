#!/bin/bash

#Partition: 'Short' (<24h) or 'Long' (>24h)
#SBATCH -p Short

#Give the job a name:
#SBATCH --job-name="summarize"

# Maximum number of tasks/CPU cores used by the job:
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5

# The amount of memory in megabytes per process in the job:
#SBATCH --mem=256000

# The maximum running time of the job in days-hours:mins:sec
#SBATCH --time=0-23:00:00

# Array
#SBATCH --array=1-38%6

# Error and output messages
#SBATCH --error=/projects/htang/CAS/C3_microbiome/2017-11-29_QIIME_CAS/log/summarize%a_%A.err
#SBATCH --output=/projects/htang/CAS/C3_microbiome/2017-11-29_QIIME_CAS/log/summarize%a_%A.out

# ----------------------
# GENERATE SUMMARY VISUALISATION OF QIIME2 SEQUENCE ARTIFACT
# ----------------------

# Set inputs
# ----------------------
src="/projects/htang/CAS/C3_microbiome/2017-11-29_QIIME_CAS/src/"
data="/projects/htang/CAS/C3_microbiome/2017-11-29_QIIME_CAS/data/"
# $1="runs_nogz.txt"
runlist=$1
line=`head -n $SLURM_ARRAY_TASK_ID $runlist | tail -1`

# Check inputs
# ----------------------

if [ "$#" -lt 1 ];
then
   echo "Missing input file location"
   echo "usage: (sbatch) script.sh <INPUT> "
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

# Run summarize
# ----------------------
# Input: SampleData[PairedEndSequencesWithQuality]
# Output: Visualisation

qiime demux summarize \
  --i-data $data/$line-demux-paired-end.qza \
  --o-visualization $data/$line-demux-paired-end.qzv

# Reset temp directory
# ----------------------
export TMPDIR=

# Exit environment
# ----------------------
source deactivate
