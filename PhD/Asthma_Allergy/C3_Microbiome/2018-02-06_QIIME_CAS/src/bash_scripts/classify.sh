#!/bin/bash

#Partition: 'Short' (<24h) or 'Long' (>24h)
#SBATCH -p Short

#Give the job a name:
#SBATCH --job-name="classify"

# Maximum number of tasks/CPU cores used by the job:
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5

# The amount of memory in megabytes per process in the job:
#SBATCH --mem=100000

# The maximum running time of the job in days-hours:mins:sec
#SBATCH --time=0-23:00:00


# Error and output messages
#SBATCH --error=/projects/htang/CAS/C3_microbiome/2017-11-29_QIIME_CAS/log/classify%J.err
#SBATCH --output=/projects/htang/CAS/C3_microbiome/2017-11-29_QIIME_CAS/log/classify%J.out

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
classifier="/projects/htang/CAS/C3_microbiome/2017-11-29_QIIME_CAS/data/gg-13-8-99-515-806-nb-classifier.qza"

# Check inputs
# ----------------------
# if [ "$#" -lt 1 ];
# then
#    echo "Missing runlist"
#    echo "usage: (sbatch) script.sh <runlist>"
#    exit 1
# fi

# Set environment
# ----------------------
source activate /usr/local/miniconda3/envs/qiime2-2017.10

# Set temp directory
# ----------------------
export TMPDIR=/projects/htang/tmp

# Classify based on pre-learnt classifier
# ----------------------
# Classify reads by taxon using a fitted classifier.

# This classifier was trained on the Greengenes 13_8 99% OTUs, where the 
# sequences have been trimmed to only include 250 bases from the region of the 
# 16S that was sequenced in this analysis (the V4 region, bound by the 
# 515F/806R primer pair). We’ll apply this classifier to our sequences, and we 
# can generate a visualization of the resulting mapping from sequence to taxonomy.

# Input: reads: FeatureData[Sequence] 
# Input: classifier: TaxonomicClassifier[Sequence] 
# Output: classification-hits: FeatureData[Taxonomy] 

qiime feature-classifier classify-sklearn \
  --i-classifier $classifier \
  --i-reads $data/rep-seqs-cas.qza \
  --o-classification $data/taxonomy-cas.qza

# Tabulate metadata (taxonomy results)
# ----------------------
# Generate a tabular view of Metadata. The output visualization supports
# interactive filtering, sorting, and exporting to common file formats.

# Input: input-file: metadata file or artifact viewable as metadata (may be multiple to merge)
# Output: Appropriate visualisation, depending on input

qiime metadata tabulate \
  --m-input-file $data/taxonomy-cas.qza \
  --o-visualization $data/taxonomy-cas.qzv

# Generate barplot from taxonomy
# ----------------------
# This visualizer produces an interactive barplot visualization of
# taxonomies. Interactive features include multi-level sorting, plot
# recoloring, category selection/highlighting, sample relabeling, and SVG
# figure export.

# Input: table: FeatureTable[Frequency]
# Input: taxonomy: FeatureData[Taxonomy] 
# Taxonomic annotations for features in the provided
# feature table. All features in the feature table
# must have a corresponding taxonomic annotation.
# Taxonomic annotations that are not present in the
# feature table will be ignored.

# Output: Appropriate visualisation

qiime taxa barplot \
  --i-table $data/table-cas.qza \
  --i-taxonomy $data/taxonomy-cas.qza \
  --m-metadata-file $data/metadata.tsv \
  --o-visualization $data/taxa-barplots-cas.qzv

# Reset temp directory
# ----------------------
export TMPDIR=

# Exit environment
# ----------------------
source deactivate
