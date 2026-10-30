#!/bin/bash

data="/Users/htang/Documents/PhD/Project - CAS/C3_Microbiome/2017-11-29_QIIME_CAS/data"

qiime tools import \
	--input-path "$data"/gg_13_5_otus/rep_set/99_otus.fasta \
	--output-path "$data"/gg_13_5_otu_99.qza \
	--type 'FeatureData[Sequence]'

qiime feature-table tabulate-seqs \
  --i-data "$data"/gg_13_5_otu_99.qza \
  --o-visualization "$data"/gg_13_5_otu_99.qzv

qiime feature-table filter-seqs \
	--i-data "$data"/rep-seqs-cas.qza \
	--i-table "$data"/table-cas.qza \
	--o-filtered-data "$data"/rep-seqs-cas-filt.qza

qiime vsearch cluster-features-closed-reference \
	--i-sequences "$data"/rep-seqs-cas-filt.qza \
	--i-table "$data"/table-cas.qza \
	--i-reference-sequences "$data"/gg_13_5_otu_99.qza \
	--p-perc-identity 1 \
	--p-threads 0 \
	--output-dir "$data"/closedref