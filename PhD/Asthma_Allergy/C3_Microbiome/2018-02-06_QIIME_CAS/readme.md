# Readme - CAS QIIME2 and Microbiome Analysis Pipeline

> **Project:** CAS (Childhood Asthma Study) microbiome analysis  
> **Data:** 16S rRNA sequencing (nasopharyngeal samples), human reference genome (for removing reads)  
> **Tools:** QIIME 2 (v2017.10), DADA2, R/phyloseq, SLURM batch scheduler  
> **Author:** Howard Tang (hhftang/hhtang)

---

## 1. Overview

This directory contains scripts for a microbiome analysis pipeline applied to the **Childhood Asthma Study (CAS)**, operated by the Telethon Kids Institute, Perth, Western Australia. Raw sequencing reads are denoised with **DADA2**, taxonomy assigned with a **classifier** trained on a bacterial reference, and downstream ecological / statistical analyses are performed in **R** using the **phyloseq** ecosystem. **QIIME2** was run in a **conda environment** on a **SLURM HPC** cluster. 

### 1.1 High-Level Flow

```
FASTQ (CASAVA 1.8 demultiplexed)
  ├──► unzip*.sh                          (Preprocessing: extract, reorganize, gzip)
  ├──► import.sh / import_simple.sh       (QIIME2: import FASTQ → .qza)
  │       ├──► dada.sh / dada_simple.sh   (QIIME2: DADA2 denoise-paired)
  │       │       └──► filthuman.sh       (QIIME2: exclude human reads)
  │       │               └──► merge_output.sh              (QIIME2: merge all samples)
  │       │                       ├──► classify.sh          (QIIME2: taxonomy + barplot)
  │       │                       ├──► summarize_output.sh  (QIIME2: summary visuals)
  │       │                       └──► export_biom.sh       (QIIME2→phyloseq export)
  │       │                               └──► R Analysis Branch [see below]
  │       ├──► summarize.sh               (QIIME2: per-sample demux viz — independent)
  │       └──► vsearch_closed_ref.sh /
  │              vsearch_closed_ref2.sh   (alternative: closed-reference OTU picking)
  └──► combine_metadata.sh                (QIIME2: combine sample metadata TSV files)

R Analysis Branch:
import_biom_format.R    (load BIOM → phyloseq)
    ├──► batch_effects.R          (DCA + PCA for batch detection)
    └──► process_phyloseq_1.R     (QC, normalization, final phyloseq objects)
            └──► core_common_microbiome.R  (core ESV/ASVs + MPG clustering)
                    ├──► compare_qiime1_qiime2_mpg.R  (QIIME1 vs QIIME2 MPG comparison)
                    ├──► age_dist.R                   (Age distribution plots)
                    └──► micro_assoc_analysis.R       (GEE: ASV/MPG vs demographics)
                            ├──► micro_assoc_analysis_vir.R   (GEE: virus associations)
                            │       ├──► time_before_and_after_infection.R  (Infection-relative labeling)
                            │       └──► phenotype_generation.R   (plink .phe files for Ober)
                            ├──► micro_global_eco.R           (Alpha diversity GEE models)
                            ├──► ssANOVA.R                    (Smooth ANOVA time series)
                            └──► hclust_MPG_trajectories.R    (Trajectory clustering)
                                    └──► rel_ab_strep.R           (Streptococcus/staphylococcus analysis)

```

---

## 2. Script index

### 2.1 Shell scripts (QIIME 2 Pipeline)

Located in `/bash_scripts/`. All scripts are designed to run on the HPC via Slurm and `sbatch`.

#### A. Data Preprocessing (Unzip) and Import

| Script | Purpose | Dependencies | Output |
| ------ | ------- | ------------ | ------ |
| `unzip[.-2].sh` | Unzip paired-end CASAVA 1.8 FASTQ, then gzip them | — | `fastq.gz` |
| `unzip_gzipped[.-3].sh`| As above, but FASTQ already gzipped | — | `fastq.gz` |
| `import.sh`| Import CASAVA 1.8 demultiplexed FASTQ → QIIME2 artifact (array: 1-38, max 6 parallel) | — | `-demux-paired-end.qza` |
| `import_simple.sh` | Simplified import (single sample, no array). | — | — |
| `summarize.sh` | Per-sample demux quality visualization (array: 1-38, max 6 parallel) | `import.sh` | `-demux-paired-end.qzv` |

#### B. Denoising and Human Read Filtering

| Script | Purpose | Dependencies | Output |
| ------ | ------- | ------------ | ------ |
| `dada.sh` | Run DADA2 denoise-paired per sample (array: 1-38, max 6 parallel). trimlength 10, trunclength 150. | `import.sh` | `-data-table.qza/qzv`, `-rep-seqs.qza/qzv` |
| `dada_simple.sh` | Simplified DADA2 (single sample, no array) | — | — |
| `filthuman.sh` |  Remove human reads via BLAST 97% identity + 97% alignment (array: 1-38, max 6 parallel) | `dada.sh` | `-filt-table.qza` |

#### C. Merge & Classification

| Script | Purpose | Dependencies | Output |
| ------ | ------- | ------------ | ------ |
| `merge_output.sh` | Merge all per-sample tables & sequences iteratively | `dada.sh`, `filthuman.sh` | `table-cas.qza`, `rep-seqs-cas.qza` |
| `summarize_output.sh` | Merged dataset summary visualization | `merge_output.sh` | `table-cas.qza`, `rep-seqs-cas.qza` |
| `classify.sh` |  Taxonomic classification using Naive Bayes classifier trained on 16S reference (Greengenes 13_8 99% OTUs) + taxa barplot | `merge_output.sh` | `taxonomy-cas.qza`, `taxa-barplots-cas.qzv` |
| `export_biom.sh` | Export QIIME2 artifacts → BIOM with taxonomy for R and phyloseq | `merge_output.sh`, `classify.sh` | `taxonomy.tsv`, `feature-table-taxa.biom` |

#### D. Other Scripts

| Script | Purpose | Dependencies | Output |
| ------ | ------- | ------------ | ------ |
| `vsearch_closed_ref[.-2].sh` | Closed-reference OTU picking with VSEARCH (separate pipeline from DADA2). | `merge_output.sh` | `clustered_sequences.qza`, `unmatched_sequences.qza`, `clustered_table.qza` |
| `combine_metadata.sh` | Combine per-sample `.tsv` metadata files into single `metadata.tsv` | — | `metadata.tsv` |

- **Note on vsearch:** DADA2 produces exact **ASVs** (amplicon sequence variants) with single-nucleotide resolution. VSEARCH closed-reference OTU picking clusters reads at a fixed % identity (typically 97%, here it is 100%) against a reference database, producing **OTUs**. These are methodologically different and produce different results. Here, the script is used to classify remaining sequences not identified with DADA2 ASVs.

---

### 2.2 R analysis scripts 

#### A. Data Import & Phyloseq Processing (QC & Normalization)

| Script | Purpose | Dependencies | Output |
| ------ | ------- | ------------ | ------ |
| `import_biom_format.R` | Load CAS metadata + QIIME2-derived BIOM feature table + taxonomy table + sample metadata CSV → phyloseq object. Sets sample IDs and taxonomy ranks. This is the **bridge** between QIIME2 shell output and all downstream R analyses. | data processing from `C0_Data_qc`, `combine_metadata.sh`, `export_biom.sh` | `import_biom_format.Rdata` |
| `batch_effects.R` | DCA + PCA ordination to check batch effects by RunName (sequencing batch), age, date, type. Uses linear models (ANOVA / linear regression) to quantify variance explained by batch and age. Generates scatter plots colored by RunName and age group.  | `import_biom_format.R` | `batch_effects.Rdata` |
| `process_phyloseq_1.R` | **Core phyloseq processing and QC:** phyloseq object → read-sum threshold qc → sum repeated samples/extractions for each log → well/mild/sick relabel and qc of logs → CSS normalization → relative abundance → core taxa extraction. | `import_biom_format.R` | `process_phyloseq_1.Rdata` |
| `get_cas_phyloseq.R` | Selects the results of a specific processing method for phyloseq, and renames the objects for later meta-analysis with COAST. | `process_phyloseq_1.R` | `cas_phyloseq_files_qc.Rdata` |

#### B. Core Microbiome & MPG Clustering

| Script | Purpose | Dependencies | Output |
| ------ | ------- | ------------ | ------ |
| `core_common_microbiome.R` | Identifies the **core microbiome** (ASV features present in ≥X% of samples) and **common microbiome** (abundance threshold >0.1%). Generates **ASV-based MPG (microbiome profile group) clusters** using **hierarchical clustering (hclust)** with Bray-Curtis distance and complete linkage. Generates barplots and heatmaps of core/common taxa . | `process_phyloseq_1.R` | `core_common_microbiome.Rdata` |
| `compare_qiime1_qiime2_mpg.R` | Code fragment; compares **MPG classifications** generated from **QIIME1** pipelines vs. **QIIME2** pipelines. | `cas_qiime1_groups.csv`, `core_common_microbiome.R`| `compare_qiime1_qiime2_mpg.Rdata` |

#### C. Association & Statistical Analyses, including (Preliminary) Trajectory Analyses

| Script | Purpose | Dependencies | Output |
| ------ | ------- | ------------ | ------ |
| `age_dist.R` | Analyzes **age distribution** of URI/LRI/Well samples. | `core_common_microbiome.R` | PDF plots |
| `micro_assoc_analysis.R` | Performs **GEE** association analyses for **general MPG types**, **season**, and **sick/well status** against ASV-level abundance (Streptococcus). Accounts for repeated measures across time. GEE models: MPG/ASV abundance vs. age, gender, visit type. | `core_common_microbiome.R` | `micro_assoc_analysis.Rdata` |
| `micro_assoc_analysis_vir.R` | Performs **GEE** association analyses specifically for **virus-level** abundance data against clinical covariates (e.g., asthma status, age, season). Accounts for repeated measures. Also explores associations/interactions with bacterial MPG. | `micro_assoc_analysis.R`  | `micro_assoc_analysis_vir.Rdata` |
| `micro_global_eco.R` | Computes **diversity metrics**: alpha diversity (Shannon, inverse Simpson) vs. age/MPG/well-sick using GEE; richness; dominance | `micro_assoc_analysis.R` | CSV outputs |
| `ssANOVA.R` | Applies **Smoothing Spline Analysis of Variance (ssANOVA)** from fitTimeSeries, fitting smooth curves to time-series microbiome abundance data (8 major ASVs over age). | `micro_assoc_analysis.R` | `ssANOVA.Rdata` |
| `hclust_MPG_trajectories.R` | Performs **hclust** of MPG trajectories over time. Precursor to MFA-K-means method seen later in `CAS_COAST_meta`. Clusters subjects by their longitudinal MPG pattern (+ silhouette optimization), then uses **GEE analysis** to test associations between cluster membership and **clinical outcomes** (URI/LRI + severity, lung function, demographics) | `micro_assoc_analysis.R` | `hclust_MPG_trajectories.Rdata` |
| `time_before_and_after_infection.R` | Analyzes microbiome composition **before vs. after infection** events, comparing healthy vs. sick states. Label healthy samples by time relative to nearest illness; visualize MPG proportion changes. Also incorporates viral data. | `micro_assoc_analysis_vir.R` | (plots) |
| `rel_ab_strep.R` | Computes **relative abundance** of key respiratory genera and ASVs: **Streptococcus, Haemophilus, Moraxella**. Tests associations (GLM) with clinical covariates (asthma status, age, season). | `hclust_MPG_trajectories.R` | — | 

#### D. Other Scripts

| Script | Purpose | Dependencies | Output |
| ------ | ------- | ------------ | ------ |
| `phenotype_generation.R` | Generates **PLINK-compatible phenotype files** from MPG classification data. Creates plink-compatible binary / categorical phenotype files for major ASVs and MPGs; for downstream genetic-wide association studies (e.g. Ober's GWAS/repeatABEL analysis) | `micro_assoc_analysis_vir.R` | `*.phe` files |

## Key Terminology

| Term | Meaning |
|------|---------|
| **ASV** | Amplicon Sequence Variant — exact biological sequences |
| **BIOM** | standardised file format for storing microbiome count and abundance data  |
| **CAS** | Childhood Asthma Study |
| **CSS** | Cumulative Sum Scaling — normalization method |
| **DADA2** | Divisive Amplicon Denoising Algorithm 2 |
| **ESV** | Exact Sequence Variant — synonym for ASV |
| **FEV1/FVC** | Forced Expiratory Volume / Forced Vital Capacity — lung function |
| **GEE** | Generalized Estimating Equations — handles repeated measures |
| **MPG** | Microbiome Profile Group — classification scheme |
| **PCoA** | Principal Coordinates Analysis |
| **PCA** | Principal Component Analysis |
| **phyloseq** | R package for microbiome analysis |
| **QIIME2** | Quantitative Insights Into Microbial Ecology (version 2) |
| **ssANOVA** | Smoothing Spline Analysis of Variance |
