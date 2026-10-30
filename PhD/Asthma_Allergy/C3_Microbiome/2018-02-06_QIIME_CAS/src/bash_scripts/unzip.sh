#!/bin/bash

#Partition: 'Short' (<24h) or 'Long' (>24h)
#SBATCH -p Short

#Give the job a name:
#SBATCH --job-name="fqunzip"

# Maximum number of tasks/CPU cores used by the job:
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5

# The amount of memory in megabytes per process in the job:
#SBATCH --mem=256000

# The maximum running time of the job in days-hours:mins:sec
#SBATCH --time=0-23:00:00

# Error and output messages
#SBATCH --error=/projects/htang/CAS/C3_microbiome/2017-11-29_QIIME_CAS/log/fqunzip_%J.err
#SBATCH --output=/projects/htang/CAS/C3_microbiome/2017-11-29_QIIME_CAS/log/fqunzip_%J.out

# ----------------------
# UNZIP CAS 16S DATA FILES 
# ----------------------

# Set inputs
# ----------------------
src="/projects/htang/CAS/C3_microbiome/2017-11-29_QIIME_CAS/src/"
data="/projects/htang/CAS/C3_microbiome/2017-11-29_QIIME_CAS/data/"
# $1="runs_nogz.txt"

# Check inputs
# ----------------------

if [ "$#" -lt 1 ];
then
   echo "Missing run list file"
   echo "usage: (sbatch) script.sh <runlist.txt>"
   exit 1
fi

# Run unzip loop
# ----------------------
while IFS='' read -r line || [[ -n "$line" ]]; do
  echo "Unzipping: $line"

# Unzip fastq
  unzip /projects/Microbiome/CAS/CAS_data/nasal/$line \
  $line/Data/Intensities/BaseCalls/*.fastq \
  -x $line/Data/Intensities/BaseCalls/flash_output/* \
  -d $data/.

	mv $data/$line/Data/Intensities/BaseCalls/* $data/$line/
	rm -r $data/$line/Data/

# Gzip fastq
	gzip $data/$line/*

# Create temporary metadata directory
  mkdir "$data"/"$line"_metadata

# Unzip samplesheet.csv
  unzip /projects/Microbiome/CAS/CAS_data/nasal/$line \
  *SampleSheet.csv \
  -d "$data"/"$line"_metadata/.

  mv "$data"/"$line"_metadata/$line/*.csv "$data"/"$line".csv
  rm -r "$data"/"$line"_metadata  

done < "$1"
