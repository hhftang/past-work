#!/bin/bash

# ----------------------
# Script to export QIIME2 outputs to something processable by phyloseq 
# ----------------------

# Problem:
# Importing qiime2-derived .biom object using read_phyloseq does not work
# due to read_phyloseq requiring taxonomy data to be incorporated in .biom file
# (which QIIME2-exported .biom format does not have)

# See:
# https://github.com/joey711/phyloseq/issues/821
# and
# https://forum.qiime2.org/t/is-there-any-way-to-summarize-taxa-plot-by-category/446/2
# for more details  

# Set inputs
# ----------------------
src="/projects/htang/CAS/C3_microbiome/2017-11-29_QIIME_CAS/src/"
data="/projects/htang/CAS/C3_microbiome/2017-11-29_QIIME_CAS/data/"

# Set environment
# ----------------------
source activate /usr/local/miniconda3/envs/qiime2-2017.10

# Export table-cas.qza to biom file
# ----------------------
qiime tools export $data/table-cas.qza \
--output-dir $data

# output file should be feature-table.biom

# Export taxonomy-cas.qza to biom file
# ----------------------
qiime tools export $data/taxonomy-cas.qza \
--output-dir $data

# output file should be taxonomy.tsv

# Edit header of taxonomy
# ----------------------

cp $data/taxonomy.tsv $data/taxonomy_0.tsv
sed 's/Feature ID/#OTUID/' $data/taxonomy_0.tsv > $data/taxonomy_1.tsv
sed 's/Taxon/taxonomy/' $data/taxonomy_1.tsv > $data/taxonomy_2.tsv
mv $data/taxonomy_2.tsv $data/taxonomy.tsv

# Add taxonomy data to .biom file
# ----------------------

biom add-metadata -i $data/feature-table.biom \
-o $data/feature-table-taxa.biom \
--observation-metadata-fp $data/taxonomy.tsv \
--sc-separated taxonomy

source deactivate