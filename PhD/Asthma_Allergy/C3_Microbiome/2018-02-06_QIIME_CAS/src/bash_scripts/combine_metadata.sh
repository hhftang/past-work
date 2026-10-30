#!/bin/bash

#Partition: 'Short' (<24h) or 'Long' (>24h)
#SBATCH -p Short

#Give the job a name:
#SBATCH --job-name="metadata.out"

# Maximum number of tasks/CPU cores used by the job:
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5

# The amount of memory in megabytes per process in the job:
#SBATCH --mem=256000

# The maximum running time of the job in days-hours:mins:sec
#SBATCH --time=0-23:00:00

# Error and output messages
#SBATCH --error=/projects/htang/CAS/C3_microbiome/2017-11-29_QIIME_CAS/log/metadata%J.err
#SBATCH --output=/projects/htang/CAS/C3_microbiome/2017-11-29_QIIME_CAS/log/metadata%J.out

# ----------------------
# CREATE COMBINED METADATA FILE
# ----------------------

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
src="/projects/htang/CAS/C3_microbiome/2017-11-29_QIIME_CAS/src/"
data="/projects/htang/CAS/C3_microbiome/2017-11-29_QIIME_CAS/data/"
runlist=$1

# Create combined metadata file
# ----------------------

cp $data/130726_M00267_0058_000000000-A48EA.tsv $data/metadata.tsv

while IFS='' read -r line; do

  sed -i -e '$a\' $data/$line.tsv
  tail -n +2 $data/$line.tsv > $data/metadata2.tsv 
	cat $data/metadata.tsv $data/metadata2.tsv > $data/metadata3.tsv
	cp $data/metadata3.tsv $data/metadata.tsv

done < $runlist

sed -i 's/Infectious/infectious/' $data/metadata3.tsv
sed -i 's/#SampleID/SampleID/' $data/metadata3.tsv

awk 'BEGIN { FS = OFS = "\t" } { for(i=1; i<=NF; i++) if($i ~ /^ *$/) $i = "NA"}; 1' $data/metadata3.tsv > $data/metadata4.tsv

sed -i $'2491s/NA\tNA\tNA\t/NA\tNA\tNA\tNA/' $data/metadata4.tsv