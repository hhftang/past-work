# CAS/COAST Microbiome Network Analysis — Notebook Inventory & Pipeline Guide

## Overview

This directory contains **26 Jupyter notebooks** (`.ipynb`) developed for the CAS and COAST microbiome cohort studies. The notebooks perform **weighted microbiome co-occurrence network construction**, **community detection** (module identification), **statistical comparison** across communities, and **distributional tests** on network topology.

---

## Quick Reference

| Category | Count | Notebooks |
|---|---|---|
| Method comparison (full datasets) | 3 | `method_test.ipynb`, `method_test_cas.ipynb`, `method_test_coast.ipynb` |
| Method comparison (subgroups) | 8 | `method_test_cas.age.leq.2.ipynb`, `method_test_cas.sick.ipynb`, `method_test_cas.well.ipynb`, `method_test_coast.age.leq.2.ipynb`, `method_test_coast.sick.ipynb`, `method_test_coast.well.ipynb` |
| Consensus clustering | 4 | `leiden_consensus.ipynb`, `leiden_RBCVP.ipynb`, `louvain_consensus.ipynb`, `Untitled.ipynb` |
| Distributional / goodness-of-fit tests | 8 | `gamma_test.ipynb`, `gamma_test_cas.ipynb`, `gamma_test_coast.ipynb`, `gamma_test_cas.age.leq.2.ipynb`, `gamma_test_cas.sick.ipynb`, `gamma_test_cas.well.ipynb`, `gamma_test_coast.age.leq.2.ipynb`, `gamma_test_coast.sick.ipynb`, `gamma_test_coast.well.ipynb` |
| Library / tool tutorials & tests | 3 | `igraph_tutorial.ipynb`, `cdlib_test.ipynb`, `cdlib.ipynb` |
| Unrelated / other | 1 | `nnt_cc.ipynb` |

> **Note:** The "subgroup" notebooks for method comparison and distributional tests cover the same cohort x subgroup matrix:
> - **CAS**: full, age ≤ 2, sick, well
> - **COAST**: full, age ≤ 2, sick, well

---

## Dependencies

All notebooks use **Python 3.8.5**.

| Library | Purpose |
|---|---|
| `igraph` | Core graph data structure; Louvain, Leiden, Infomap, CPM, SLM community detection |
| `leidenalg` | Leiden algorithm implementation (with multi-level refinement) |
| `cdlib` | Unified community detection library wrapper (LFM, OSLOM, greedy modularity, CPM, etc.) |
| `networkx` | Alternative graph library; some conversion / utility functions |
| `numpy` | Numerical operations, random sampling |
| `matplotlib` | Plotting (histograms, heatmaps, bar charts, line plots) |
| `seaborn` | Statistical visualisation, heatmaps |
| `cairo` | Graphical rendering for igraph |
| `palettable` | Color palettes |
| `scipy.stats` | Statistical tests (gamma distribution fitting, KS tests, chi-square) |
| `netneurotools` | **Only in `nnt_cc.ipynb`** — neuroimaging network analysis |

---

## Inputs & Outputs

### Primary Input File

```
/Users/htang/Documents/PhD/Project - CAS/C3_Microbiome/2019-10-24_QIIME_CAS_COAST_meta_5/output/2020-09-25_fastspar_community/net.graph.all.pos.txt
```

- **Format:** GML (Graph Markup Language)
- **Nodes:** 532 (OTUs / microbial taxa)
- **Edges:** 23,274 (significant positive SparCC correlations)
- **Cohort subsetting** is performed *in code* by filtering nodes whose metadata columns (age, health status) match the desired subgroup.

### Outputs (generated files)

All notebooks write output to the same base directory under `/Users/htang/…`:

```
/Users/htang/Documents/PhD/Project - CAS/C3_Microbiome/2019-10-24_QIIME_CAS_COAST_meta_5/output/2020-09-25_fastspar_community/
```

| Output | Description |
|---|---|
| `net.graph.all.pos.txt` | Input GML network |
| `method_test_cas.txt` / `method_test_coast.txt` | Community membership assignments per method |
| `consensus_heatmap_*_cas.png` / `*_coast.png` | Consensus co-assignment heatmaps |
| `membership_plot_*_cas.png` / `*_coast.png` | Membership distribution bar charts |
| `gamma_fit_cas.png` / `*_coast.png` | Gamma distribution fit plots on node degrees |
| Various `fig` PNGs | Plots from igraph tutorial and cdlib tests |

---

## Notebook-by-Notebook Summary

### 1. `igraph_tutorial.ipynb`
- **Purpose:** Learning / demonstration notebook for `igraph` — basic graph operations, visualization, and introduction to the Leiden algorithm.
- **Input:** GML network file.
- **Output:** Plotting examples; no analytical outputs.
- **Role in pipeline:** Foundational reference; not used in any analysis downstream.

### 2. `cdlib_test.ipynb`
- **Purpose:** Exploratory testing of `cdlib` library functions — verifying which community detection algorithms are available and how they behave.
- **Input:** GML network file.
- **Output:** Diagnostic plots; no formal results.
- **Role in pipeline:** Exploratory / validation; standalone.

### 3. `cdlib.ipynb`
- **Purpose:** Applying `cdlib` wrapper to run all supported community detection methods on the full CAS and COAST datasets. Outputs community labels per node.
- **Input:** GML network file.
- **Output:** `method_test_cas.txt` / `method_test_coast.txt` (or equivalent label files).
- **Role in pipeline:** Core analysis — replaces/parallels the `method_test*.ipynb` family.

### 4. `method_test.ipynb`
- **Purpose:** Method comparison on the **full combined dataset**. Runs Louvain, Leiden, Infomap, LFM, SLM, OSLOM, greedy modularity, and CPM. Generates community membership plots and consensus heatmaps.
- **Input:** GML network file.
- **Output:** Membership plots and consensus heatmaps for the full network.
- **Role in pipeline:** Baseline full-dataset comparison.

### 5. `method_test_cas.ipynb`
- **Purpose:** Identical pipeline to `method_test.ipynb` but restricted to the **CAS cohort** only.
- **Input:** GML network file (subsetted to CAS metadata).
- **Output:** Method comparison outputs for CAS.
- **Role in pipeline:** Primary cohort comparison.

### 6. `method_test_coast.ipynb`
- **Purpose:** Identical pipeline to `method_test.ipynb` but restricted to the **COAST cohort** only.
- **Input:** GML network file (subsetted to COAST metadata).
- **Output:** Method comparison outputs for COAST.
- **Role in pipeline:** Secondary cohort comparison.

### 7. `method_test_cas.age.leq.2.ipynb`
- **Purpose:** Method comparison on **CAS children aged ≤ 2 years**.
- **Input:** GML network file (subsetted to CAS, age ≤ 2).
- **Output:** Method comparison outputs for CAS age ≤ 2.
- **Role in pipeline:** Subgroup analysis.

### 8. `method_test_cas.sick.ipynb`
- **Purpose:** Method comparison on **CAS subjects with sickness/illness**.
- **Input:** GML network file (subsetted to CAS, sick).
- **Output:** Method comparison outputs for CAS sick.
- **Role in pipeline:** Subgroup analysis.

### 9. `method_test_cas.well.ipynb`
- **Purpose:** Method comparison on **CAS healthy/well subjects**.
- **Input:** GML network file (subsetted to CAS, well).
- **Output:** Method comparison outputs for CAS well.
- **Role in pipeline:** Subgroup analysis.

### 10. `method_test_coast.age.leq.2.ipynb`
- **Purpose:** Method comparison on **COAST children aged ≤ 2 years**.
- **Input:** GML network file (subsetted to COAST, age ≤ 2).
- **Output:** Method comparison outputs for COAST age ≤ 2.
- **Role in pipeline:** Subgroup analysis.

### 11. `method_test_coast.sick.ipynb`
- **Purpose:** Method comparison on **COAST sick subjects**.
- **Input:** GML network file (subsetted to COAST, sick).
- **Output:** Method comparison outputs for COAST sick.
- **Role in pipeline:** Subgroup analysis.

### 12. `method_test_coast.well.ipynb`
- **Purpose:** Method comparison on **COAST healthy/well subjects**.
- **Input:** GML network file (subsetted to COAST, well).
- **Output:** Method comparison outputs for COAST well.
- **Role in pipeline:** Subgroup analysis.

### 13. `leiden_consensus.ipynb`
- **Purpose:** Runs **Leiden algorithm 100 times** with different random seeds on full CAS and COAST datasets. Computes majority-vote **consensus membership** (how often two nodes are co-clustered). Produces consensus heatmaps and membership plots.
- **Input:** GML network file.
- **Output:** `consensus_heatmap_*_cas.png`, `membership_plot_*_cas.png` (and COAST equivalents).
- **Role in pipeline:** Robustness assessment — Leiden stability analysis.

### 14. `leiden_RBCVP.ipynb`
- **Purpose:** Similar 100-iteration Leiden consensus pipeline, but with **custom plotting functions** specific to the RBCVP (Research Biobank Core/Vascular Phenotyping) project. May include different visualization styles or additional annotations.
- **Input:** GML network file.
- **Output:** Custom Leiden consensus visualizations.
- **Role in pipeline:** Robustness assessment — Leiden stability, project-specific visualization.

### 15. `louvain_consensus.ipynb`
- **Purpose:** Runs **Louvain algorithm 100 times** with different random seeds. Computes consensus membership and produces heatmaps / membership plots. Allows direct comparison of Louvain vs. Leiden stability.
- **Input:** GML network file.
- **Output:** `consensus_heatmap_louvain_*_cas.png`, `membership_plot_louvain_*_cas.png` (and COAST equivalents).
- **Role in pipeline:** Robustness assessment — Louvain stability analysis.

### 16. `Untitled.ipynb`
- **Purpose:** Leiden 100-iteration consensus analysis with **custom plotting functions** defined inline. Appears to be a variant/branch of the consensus notebooks with bespoke visualization logic.
- **Input:** GML network file.
- **Output:** Custom consensus visualizations.
- **Role in pipeline:** Robustness assessment — Leiden stability, ad-hoc plotting.

### 17. `gamma_test.ipynb`
- **Purpose:** Tests whether the **node degree distribution** of the full combined network follows a **gamma distribution**. Fits gamma parameters and produces diagnostic plots (histogram + fitted curve).
- **Input:** GML network file.
- **Output:** `gamma_fit.png` — goodness-of-fit plot.
- **Role in pipeline:** Topological characterization — standalone distribution test.

### 18. `gamma_test_cas.ipynb`
- **Purpose:** Same gamma distribution test on the **CAS cohort** subnetwork.
- **Input:** GML network file (subsetted to CAS).
- **Output:** `gamma_fit_cas.png`.
- **Role in pipeline:** Topological characterization.

### 19. `gamma_test_coast.ipynb`
- **Purpose:** Same gamma distribution test on the **COAST cohort** subnetwork.
- **Input:** GML network file (subsetted to COAST).
- **Output:** `gamma_fit_coast.png`.
- **Role in pipeline:** Topological characterization.

### 20. `gamma_test_cas.age.leq.2.ipynb`
- **Purpose:** Gamma distribution test on **CAS age ≤ 2** subnetwork.
- **Input:** GML network file (subsetted to CAS, age ≤ 2).
- **Output:** `gamma_fit_cas.age.leq.2.png` (or equivalent).
- **Role in pipeline:** Subgroup topological characterization.

### 21. `gamma_test_cas.sick.ipynb`
- **Purpose:** Gamma distribution test on **CAS sick** subnetwork.
- **Input:** GML network file (subsetted to CAS, sick).
- **Output:** `gamma_fit_cas.sick.png`.
- **Role in pipeline:** Subgroup topological characterization.

### 22. `gamma_test_cas.well.ipynb`
- **Purpose:** Gamma distribution test on **CAS well** subnetwork.
- **Input:** GML network file (subsetted to CAS, well).
- **Output:** `gamma_fit_cas.well.png`.
- **Role in pipeline:** Subgroup topological characterization.

### 23. `gamma_test_coast.age.leq.2.ipynb`
- **Purpose:** Gamma distribution test on **COAST age ≤ 2** subnetwork.
- **Input:** GML network file (subsetted to COAST, age ≤ 2).
- **Output:** `gamma_fit_coast.age.leq.2.png`.
- **Role in pipeline:** Subgroup topological characterization.

### 24. `gamma_test_coast.sick.ipynb`
- **Purpose:** Gamma distribution test on **COAST sick** subnetwork.
- **Input:** GML network file (subsetted to COAST, sick).
- **Output:** `gamma_fit_coast.sick.png`.
- **Role in pipeline:** Subgroup topological characterization.

### 25. `gamma_test_coast.well.ipynb`
- **Purpose:** Gamma distribution test on **COAST well** subnetwork.
- **Input:** GML network file (subsetted to COAST, well).
- **Output:** `gamma_fit_coast.well.png`.
- **Role in pipeline:** Subgroup topological characterization.

### 26. `nnt_cc.ipynb` ⚠️
- **Purpose:** **Unrelated to the microbiome pipeline.** Downloads structural brain connectivity data from OpenNeuro (dataset ds000031) and analyzes it using `netneurotools` for neuroimaging connectomics research. Extracts cortical thickness (CT), sulcal depth (SD), and their covariance across brain regions.
- **Input:** Downloads from `https://openneuro.org/datasets/ds000031`.
- **Output:** Neuroimaging brain network analyses.
- **Role in pipeline:** **None** — separate research project stored in the same directory.

---

## Data Sensitivity & Personal Information

| Item | Location | Recommendation |
|---|---|---|
| **Hardcoded macOS path** | `/Users/htang/Documents/PhD/Project - CAS/…` | Found in ~20 notebooks. Replace with a relative path, environment variable, or a path resolved at runtime (e.g., `os.path.expanduser("~")` or a config file). |
| **macOS username** | `htang` | Embedded in all file paths above. Consider parameterising. |
| **Cohort metadata columns** | Used for subsetting (age, health/sickness status) | The raw metadata CSV/TSV file is referenced implicitly by the path above. Ensure it is not committed to version control if it contains PII. |
| **No passwords/passcodes found** | — | — |

---

## Inferred Pipeline Execution Order

```
Phase 0 — Exploration & Tooling
├── igraph_tutorial.ipynb          (learn igraph basics)
├── cdlib_test.ipynb               (test cdlib API)
├── cdlib.ipynb                    (run cdlib on full datasets)
│
Phase 1 — Method Comparison (community detection)
│   method_test.ipynb              (full combined network)
│   ├── method_test_cas.ipynb      (CAS cohort)
│   │   ├── method_test_cas.age.leq.2.ipynb
│   │   ├── method_test_cas.sick.ipynb
│   │   └── method_test_cas.well.ipynb
│   └── method_test_coast.ipynb    (COAST cohort)
│       ├── method_test_coast.age.leq.2.ipynb
│       ├── method_test_coast.sick.ipynb
│       └── method_test_coast.well.ipynb
│
Phase 2 — Consensus / Robustness (100-iteration clustering)
│   leiden_consensus.ipynb         (Leiden, full datasets)
│   └── leiden_RBCVP.ipynb         (Leiden, project-specific plots)
│   louvain_consensus.ipynb        (Louvain, full datasets)
│   Untitled.ipynb                 (Leiden, ad-hoc plots)
│
Phase 3 — Topological Characterisation (degree distribution)
│   gamma_test.ipynb               (full combined network)
│   ├── gamma_test_cas.ipynb
│   ├── gamma_test_coast.ipynb
│   │   ├── gamma_test_cas.age.leq.2.ipynb
│   │   ├── gamma_test_cas.sick.ipynb
│   │   ├── gamma_test_cas.well.ipynb
│   │   ├── gamma_test_coast.age.leq.2.ipynb
│   │   ├── gamma_test_coast.sick.ipynb
│   │   └── gamma_test_coast.well.ipynb
│
Phase X — Unrelated
└── nnt_cc.ipynb                   (neuroimaging, not part of this pipeline)
```

### Key Observations

1. **Phase 1 and Phase 3 are parallel and independent** — they both subset from the same source network but do not depend on each other's outputs.
2. **Phase 2 (consensus) follows Phase 1** in spirit: once you know which algorithms are being used, you assess their stability. However, there is no hard dependency — each consensus notebook is self-contained.
3. **`cdlib.ipynb` may be an earlier version** of `method_test*.ipynb` family; both produce similar outputs (community labels). Consider consolidating.
4. **Four consensus notebooks** (`leiden_consensus`, `leiden_RBCVP`, `louvain_consensus`, `Untitled`) produce overlapping outputs. The main difference is the algorithm (Leiden vs. Louvain) and the visualization style. These could be merged into a single configurable notebook.
5. **`nnt_cc.ipynb` has no relationship** to the microbiome analysis and should ideally be moved to a separate directory.

---

## Suggested Consolidation

To reduce redundancy, consider merging the notebook family into a **single parameterised pipeline**:

```python
# Pseudocode for a consolidated notebook
cohorts = ["cas", "coast"]
subgroups = ["full", "age_leq_2", "sick", "well"]
methods = ["louvain", "leiden", "infomap", "lfm", "slm", "oslom", "greedy", "cpm"]
iterations = [1, 100]  # single-run vs. consensus

for cohort in cohorts:
    for subgroup in subgroups:
        G = load_network(cohort, subgroup)
        for method in methods:
            labels = run_community_detection(G, method)
            save_labels(labels, f"{cohort}_{subgroup}_{method}.txt")
        for n_iter in [100]:
            labels_list = [run_community_detection(G, "leiden", seed=i) for i in range(100)]
            consensus = compute_consensus(labels_list)
            plot_consensus_heatmap(consensus, f"consensus_{cohort}_{subgroup}.png")
        degree = G.degree()
        plot_gamma_fit(degree, f"gamma_{cohort}_{subgroup}.png")
```

---

*Generated from analysis of all 26 `.ipynb` files in this directory.*
