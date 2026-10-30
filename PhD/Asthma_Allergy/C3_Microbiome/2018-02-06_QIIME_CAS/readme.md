# Readme - QIIME2 in CAS

## Overview

This directory contains scripts for a microbiome analysis pipeline applied to the **Childhood Asthma Study (CAS)**, Perth, Western Australia. Raw sequencing reads are denoised with **DADA2**, taxonomy assigned with **QIIME2**, and downstream ecological / statistical analyses are performed in **R** using the **phyloseq** ecosystem.

Uses `qiime2-2017.10` **conda environment** on a SLURM HPC cluster.

Requires **human reference genome** (for removing host reads).

## Bash scripts

Located in `/bash_scripts/`. All scripts are designed to run on the HPC via Slurm and `sbatch`.

In order:

1. `unzip[.-3].sh`: Unzip paired-end `.fastq.gz` files (R1/R2) from a manifest list.
2. `unzip_gzipped[.-3].sh`: Unzip successive batches of gzipped FASTQ files (single-end or combined)
3. `import_simple.sh`: Import paired FASTQ files into a QIIME2 artifacts directory using `qiime tools import` with `PairedEndFASTQV2`.
    - `import.sh`: Import paired FASTQ files (alternate version; including manifest-based import).
4. `dada_simple.sh`: Run DADA2 denoising (`qiime dada2 denoise-paired`) with simplified parameters and fewer QC checks. 
    - `dada.sh`: Run DADA2 denoising (`qiime dada2 denoise-paired`) with full/standard parameters and multiple QC steps.
5. `summarize_output.sh`: Summarize DADA2 output (sequence-length plots, merge stats, etc.). 
6. `merge_output.sh`: Merge forward and reverse reads post-DADA2 (`qiime feature-table merge`). 
7. `classify.sh`: Assign taxonomy using `qiime feature-classifier classify-sklearn` (Naive Bayes classifier on 16S reference). 
8. `filthuman.sh`: Filter out human (host) contaminating reads using VSEARCH (i.e. `qiime feature-classifier classify-consensus-vsearch`) matched against human reference genome.
9. `combine_metadata.sh`: Merge sample metadata (phenotype / clinical data) with QIIME2 feature table / taxonomy output.
10. `summarize.sh`: Generate QIIME2 summary visualizations (e.g., taxa bar plots, phylogenetic trees).
11. `export_biom.sh`: Export QIIME2 artifacts (feature table, taxonomy) to **BIOM** format for downstream R import.

Other scripts:

- `vsearch_closed_ref[.-2].sh`: Closed-reference OTU picking with VSEARCH (separate pipeline from DADA2). DADA2 produces exact **ASVs** (amplicon sequence variants) with single-nucleotide resolution. VSEARCH closed-reference OTU picking clusters reads at a fixed % identity (typically 97%) against a reference database, producing **OTUs**. These are methodologically different and produce different results.

- `dada_sum.sh`: Summary-only: reads DADA2 output files and produces a summary report (assumes DADA2 already ran). 

## R scripts 

Core processing scripts - In order:

1. `import_biom_format.R`: Loads a QIIME2-derived BIOM feature table + taxonomy table + sample metadata CSV; converts to a `phyloseq` object. Sets sample IDs and taxonomy ranks. This is the **bridge** between QIIME2 shell output and all downstream R analyses.

2. `process_phyloseq_1.R`: Core phyloseq processing: filters rare features, rarefies the feature table, computes **CSS** (Cumulative Sum Scaling) normalization, calculates **alpha diversity** (Shannon, Observed, Faith's PD) and **beta diversity** (Bray-Curtis, UniFrac), performs **PCoA** and **PCA** ordination, and generates diversity plots. Saves the processed `phyloseq` object as `.Rdata`.

3. (opt.) `get_cas_phyloseq.R`: Loads the intermediate `phyloseq` object from `2018-02-06_process_phyloseq_1.Rdata`, applies a final QC step (filtering out samples < 3000 reads), and saves the final `cas_phyloseq_files_qc.Rdata`. This is the **last QC step** before all downstream analyses.

Downstream analyses:

- `batch_effects.R`: Tests whether sequencing batch (**RunName**) and **age** affect the first two PCoA axes (PC1, PC2). Uses linear models (ANOVA / linear regression) to quantify variance explained by batch and age. Generates scatter plots colored by RunName and age group.

- `micro_global_eco.R`: Computes **global ecosystem diversity metrics**: alpha diversity (Shannon, Simpson, Observed), beta diversity (Bray-Curtis PCoA), and community dominance/evenness across all samples. 

- `ssANOVA.R`:  Applies **Smoothing Spline Analysis of Variance (ssANOVA)** curve fitting to time-series microbiome abundance data. Fits smooth curves to genus-level abundance over time, identifies significant temporal trends.

- `micro_assoc_analysis_vir.R`: Performs **GEE (Generalized Estimating Equations)** association analyses specifically for **virus-level** abundance data against clinical covariates (e.g., asthma status, age, season). Accounts for repeated measures.

- `core_common_microbiome.R`: Identifies the **core microbiome** (features present in ≥X% of samples) and **common microbiome** (abundance threshold). Generates **ESV-based MPG (microbiome profile group) classification**. Generates bar plots of core/common taxa.

- `phenotype_generation.R`: Generates **PLINK-compatible phenotype files** from MPG classification data. Creates binary / categorical phenotype files for downstream genetic / microbiome-wide association studies.

- `time_before_and_after_infection.R`: Analyzes microbiome composition **before vs. after infection** events, comparing healthy vs. sick states. Uses paired statistical tests (Wilcoxon / t-test) and PCoA visualization.

- `age_dist.R`: Analyzes **age distribution** of samples across cohorts, MPG groups, and clinical categories. Generates histograms, density plots, and statistical comparisons.

- `micro_assoc_analysis.R`: Performs **GEE (Generalized Estimating Equations)** association analyses for **general MPG types**, **season**, and **sick/well status** against genus-level and ESV-level abundance. Accounts for repeated measures across time.

- `hclust_MPG_trajectories.R`: Performs **hierarchical clustering** of MPG trajectories over time. Clusters subjects by their longitudinal MPG pattern, then uses **GEE analysis** to test associations between cluster membership and **clinical outcomes** (FEV1, FVC, FEV1/FVC ratio).

- `rel_ab_strep.R`: Computes **relative abundance** of key respiratory genera: **Streptococcus, Haemophilus, Moraxella**. Tests associations with clinical covariates (asthma status, age, season). Generates bar plots and volcano plots.

- `compare_qiime1_qiime2_mpg.R`: Compares **MPG classifications** generated from **QIIME1** pipelines vs. **QIIME2** pipelines. Computes concordance metrics (Cohen's kappa, confusion matrix). Generates scatter plots and agreement tables.

## Key Terminology

| Term | Meaning |
|------|---------|
| **CAS** | Childhood Asthma Study |
| **DADA2** | Divisive Amplicon Denoising Algorithm 2 |
| **ASV** | Amplicon Sequence Variant — exact biological sequences |
| **ESV** | Exact Sequence Variant — synonym for ASV |
| **QIIME2** | Quantitative Insights Into Microbial Ecology (version 2) |
| **BIOM** | standardised file format for storing microbiome count and abundance data  |
| **phyloseq** | R package for microbiome analysis |
| **CSS** | Cumulative Sum Scaling — normalization method |
| **MPG** | Microbiome Profile Group — classification scheme |
| **GEE** | Generalized Estimating Equations — handles repeated measures |
| **ssANOVA** | Smoothing Spline Analysis of Variance |
| **PCoA** | Principal Coordinates Analysis |
| **PCA** | Principal Component Analysis |
| **FEV1/FVC** | Forced Expiratory Volume / Forced Vital Capacity — lung function |
