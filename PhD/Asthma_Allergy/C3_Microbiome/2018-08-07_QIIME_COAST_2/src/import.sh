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

# Array
#SBATCH --array=1-41%6

# Error and output messages
#SBATCH --error=/projects/htang/COAST/2017-12-08_QIIME_COAST/log/fqimport%a_%A.err
#SBATCH --output=/projects/htang/COAST/2017-12-08_QIIME_COAST/log/fqimport%a_%A.out

# ----------------------
# IMPORT FASTQ FILES TO QIIME2 ARTIFACT FORMAT
# ----------------------

# Set inputs
# ----------------------
src="/projects/htang/COAST/2017-12-08_QIIME_COAST/src/"
data="/projects/htang/COAST/2017-12-08_QIIME_COAST/data/"
runlist=$1
line=`head -n $SLURM_ARRAY_TASK_ID $runlist | tail -1`

# Check inputs
# ----------------------

if [ "$#" -lt 1 ];
then
   echo "Missing run list file"
   echo "usage: (sbatch) script.sh <runlist.txt>"
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
