# Survival-XML-TME-BC

This repository contains the code and data supporting the manuscript *Integrated multiomics analysis unveils how macrophages drive immune suppression in breast tumors and affect clinical outcomes*. This study investigates the role of immune-suppressive macrophages in breast tumors, integrating multiomics data such as bulk RNA-seq, scRNA-seq and imaging mass cytometry (IMC) data. Follow the steps below to recreate the study’s results.

## Overview of Analysis Steps

### 1. Estimating Cell Fractions
To estimate cell fractions for bulk RNA-seq profiles, use the Signature Matrices available in our previous work: https://github.com/YounessAzimzade/XML-TME-NAC-BC. We recommend using **CIBERSORTx** for cell fraction estimation. After estimating fractions and retrieving clinical outcome data, you should have datasets similar to `NAC.csv` for the NAC cohort and `MBRC.csv`/`TCGA.csv` for other cohorts.

### 2. Analyzing the Role of Cell Types in Clinical Outcomes
Explore the influence of cell type frequencies on clinical outcomes using **Survival Support Vector Machines (SSVM)** and **Random Survival Forests (RSF)**:
   - **SSVM.ipynb** and **RSF.ipynb**: Train models on cell fractions, clinical features, and outcomes to predict clinical outcomes.
   - **SHAP analysis**: Extract feature importance using SHAP values, yielding SHAP datasets with processed feature columns (performed in SSVM.ipynb and RSF.ipynb).  
   - After SHAP analysis, calculate "Survival Scores" using **fig2.R**.

### 3. Juxtaposing Survival and pCR Scores
Once "Survival Scores" are calculated, these scores are compared with pathological complete response (pCR) scores (available at https://github.com/YounessAzimzade/XML-TME-NAC-BC) using **fig3.R**.

### 4. Classical Survival Analysis
For additional insights, perform traditional survival analyses using scripts in the **Survival Analysis** folder.

### 5. Imaging Mass Cytometry (IMC) and scRNA-seq Analysis
Analyze IMC and single-cell RNA-seq data using scripts in their respective folders. Follow folder-specific instructions to reproduce analyses and generate visualizations for spatial organization and immune profiling insights.

## Repository Structure

- **data/**  
  Includes files for cell fractions calculated using deconvolution (TCGA and MBRC), cell fractions from IMC data, distances calculated for macrophages and other cells in IMC samples as well as macrophages and T cells from scRNA-seq data and cell frequencies in scRNA-seq samples.

- **scripts/**  
  Contains scripts for each stage of analysis:
  - **SSVM.ipynb** and **RSF.ipynb**: Train models to predict clinical outcomes and compute SHAP values.
  - **fig2.R**: Calculates Survival Scores from SHAP values.
  - **fig3.R**: Compares Survival Scores with pCR scores.
  - **IMC**: Spatial analysis of cell phenotypes with cell-cell distance calculations.
    - Calculate min distance macs vs all.R
    - Cell Type Fractions.R
    - Fig 4f.R
    - Fig 4g.R
    - Fig 4h.R
    - Fig 5a.R
    - Fig 5c.R
    - Fig 5d & e.R
  - **scRNA-seq**: Analyzes cell-frequency correlations and macrophage-T cell interactions using NicheNet.
    - Annotation Transfer.R
    - Fig5 g & h.R
    - Fig5 h, i, j & k.R
    - Fig5 l, m, n, o.R

## Data Access

Data sources:
- **TCGA** and **MBRC** datasets: Available on cBioPortal (https://www.cbioportal.org/).
- **Spatial omics data**:
  - Danenberg et al., 2022: https://zenodo.org/record/5850952
  - Wang et al., 2023: https://zenodo.org/records/7990870
- **scRNA-seq data**: 
  - Wu et al., 2021: https://singlecell.broadinstitute.org/single_cell/study/SCP1039/a-single-cell-and-spatially-resolved-atlas-of-human-breast-cancers
  - Bassez et al., 2021: https://lambrechtslab.sites.vib.be/en/single-cell

## Contact

For questions or support, please contact Youness Azimzade at younessazimzade@gmail.com.
