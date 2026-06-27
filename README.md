# Smad1 Protects the Remodeling Infarcted Heart, Stimulating a Fibrosis-Restraining Myofibroblast Subpopulation. Hanna et al., Circulation, 2026.
**scRNA-seq Integrative Analysis Pipeline**

This repository contains the code required to reproduce the single-cell RNA sequencing (scRNA-seq) integration, clustering, and functional enrichment analysis for the manuscript: 
*"Smad1 Protects the Remodeling Infarcted Heart, Stimulating a Fibrosis-Restraining Myofibroblast Subpopulation"* (Hanna et al., Circulation, 2026).

## Overview
To investigate cardiac fibroblast heterogeneity and Smad1 pathway activity during myocardial infarction healing, we performed an integrative analysis of two publicly available scRNA-seq datasets. 

## Requirements & Environment
This pipeline was built and executed in R. The integration specifically utilizes **Seurat v5** native objects. A complete list of package versions and system specs can be found in the `session_info.txt` file.

**Core Dependencies:**
* Seurat (v5.0+)
* harmony
* SCP
* clusterProfiler

## Usage
1. Clone this repository to your local machine.
2. Ensure the raw matrices (from GEO/ArrayExpress) are in the project folder.
3. Run `Hanna_Circulation_2026.R` sequentially.
