# Readme - C0 - Data quality control

## 2017-08-04

Contains attempts to clean CAS datasets and harmonise phenotypes so that they
are (somewhat) internally-consistent with each other.

The most troublesome aspect of data cleaning and quality control concerns
the viral sampling and respiratory infection phenotypes. This is because of:

- Differences in viral typing methods (Johnson vs. Gern Lab for Year 1 samples)
- Differences in definition of respiratory infection phenotypes (e.g. 
between CAS and COAST) 

There will be a meeting at some stage in the near feature to try to get
a set of consensus phenotypes out, that *everyone* (Inouye Lab and Telethon Labs) 
can use. 

# 2026-10-30

## Main scripts

In order:

1. `immuno_var_gen.R`: Loads immunology datasets for ages 0–5 and age 10, then merges them. Loads CMI (cell-mediated immunity) simulation data from seven time points (y0, y0.5, y1, y2, y3, y5) and merges. Loads Poly-IC stimulation data. Merges everything into a single dataframe immunol_cyto_all.df. Performs quality checks comparing dob/DOB.x/DOB.y (date of birth) and sex/GENDER.x/GENDER.y (gender) across datasets, then drops the redundant columns.

2. `vir_respchart_var_gen.R`: Loads viral pathogen data (NPA samples for years 1, 2, 3) and HRV subtyping data. Loads clinical ARI data (CAS_PNA_cyrobanked_171013.txt). Merges these per year. Loads the respchart dataset and generates new phenotype variables (uri_new, cari, fari). Merges viral and respchart datasets, performs consistency checks, generates per-individual summary counts at monthly (0–2m, 2–6m, 6–12m) and yearly (1–5) timepoints, respecting exclusion ages from the immunology dataset. Computes summary statistics (sum, mean, median) and generates PDF histograms for URI, LRI, fari, cari, and wlri phenotypes. Additionally adds new phenotypes: sari (severe ARI = fever OR wheeze), slri (severe LRI = (fever OR wheeze) AND LRI), and wflri (wheezy+febrile LRI), and creates new timepoint 425d (year 1 = 0–425 days). Generates histograms for all new phenotype types. Finally merges with immunology data and performs extensive cross-dataset consistency checks between the new vir/respchart counts and the older immunology-derived counts, including detailed per-subject breakdowns of HRV and RSV discrepancies with multiple data frame extractions. 

3. `virus_baseline_timeframe_check.R`: Checks whether any baseline (control) viral samples were collected within 28 days (4 weeks) after an infection sample for the same subject — i.e., potential contamination of baseline measurements by a recent infection. Performs this check separately for year 1, year 2, year 3, and the full clinical ARI dataset. Outputs counts of problematic baseline–infection pairs.

4. `excl_age_doublecheck.R`: Re-computes infection counts from virus_resp.df without applying the exclusion age filter, producing versions of the summary dataframes with noexcl_ in the name. Then applies exclusion age filters and compares the two sets. Also reads an external CSV file (data/new ari variables.csv) containing new exclusion age data and merges it. Designed to double-check that the exclusion logic in the earlier scripts is correct.

5. `master_table_gen.R`: Loads the vir/respchart output from `2016-10-26_vir_respchart_var_gen.Rdata`. Renames columns (e.g. hdmigg10 → hdmigg_10). Selects a curated subset of variables spanning infection characteristics, allergen IgE/IgG/IgG4, SPT wheal sizes, family history, environmental exposures, flow cytometry counts, TREC counts, body size, breast milk composition, maternal vitamin D, mRNA/protein expression data, and clinical phenotypes. Applies log10 transformations to antibody assays, TRECS, mRNA, and protein values (censoring zeros to small positive values first). Collapses cat exposure from three levels (none/indoor/outdoor) to binary (none/yes). Writes column lists to text files and saves both a full master table (with SUBJ_ID) and a numeric-only master table.

## Other scripts

- `special_wf_tables.Rmd`: RMarkdown document that: (1) defines wheeze/fever ARI phenotypes (lri, wlri, fari, flri, wflri, sari, slri) and timepoints; (2) loads the Rdata from a prior vir/respchart run and displays summary tables of total counts and mean counts per individual; (3) generates PDF histograms for LRI, WLRI, FARI, FLRI, WFLRI, SARI, and SLRI phenotypes across age groups.

- `get_slri.R`: Single one-liner: extracts SUBJ_ID, slri_1 through slri_5 from premaster.df and writes to output/slri.txt.

- `y10_data_present.R`: Generates a series of heatmaps showing data completeness (presence/absence of non-NA values) across subjects for year-10 data from three datasets: (1) immunology (immunol_10.df and related), (2) polyIC stimulation (polyic_y4_10.df), and (3) lung function/SPT (10yr LFT and SPT.csv). Uses the GMD package's heatmap.3 with a warm colour palette. Produces a final composite heatmap and writes results to CSV.

**NB: Log-transform caveats:** The scripts note that censoring zero values (to 0.000001 for mRNA, 0.01 for protein) before log10 transformation is suboptimal, and obtaining raw data is recommended. See lines 438–447 of the master_table_gen scripts and line 40 in the immuno_var_gen script.