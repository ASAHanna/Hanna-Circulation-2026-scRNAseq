# Smad1 Protects the Remodeling Infarcted Heart, Stimulating a Fibrosis-Restraining Myofibroblast Subpopulation. Hanna et al.
**scRNA-seq Integrative Analysis Pipeline**

This repository contains the code required to reproduce the single-cell RNA sequencing (scRNA-seq) integration, clustering, and functional enrichment analysis for the manuscript: 
*"Smad1 Protects the Remodeling Infarcted Heart, Stimulating a Fibrosis-Restraining Myofibroblast Subpopulation"*

## Overview
To comprehensively investigate cardiac fibroblast heterogeneity and Smad1 pathway activity during myocardial infarction healing, we performed an integrative analysis of two distinct publicly available scRNA-seq datasets. 

**Methodological Advantage & Statistical Power:**
By integrating data from two distinct fibroblast reporter models (*Pdgfra*-EGFP and *Col1a1*-GFP), this pipeline drastically increases the cellular depth at the critical 7-day post-MI timepoint—the peak of the proliferative phase with the highest myofibroblast density. 

Because the clustering space utilizes batch-corrected parameters derived from this combined transcriptomic pool, the global unsupervised clustering naturally possesses the statistical power to resolve the highly activated myofibroblast compartment into distinct functional states (Regulatory Myofibroblasts [RMFs] and Fibrogenic Myofibroblasts [FMFs]) directly on the primary manifold, bypassing the need for secondary isolation or recursive sub-clustering. 

## Requirements & Environment
This pipeline was built and executed in R. A complete list of package versions and system specs can be found in the `session_info.txt` file.

## Usage
1. Clone this repository to your local machine.
2. Ensure the raw matrices (from GEO/ArrayExpress) are in the project folder.
3. Run `Hanna_Circulation_2026.R` sequentially.
