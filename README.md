# Survival-XML-TME-BC

This repository contains the code and data supporting the manuscript *Integrated Multiomics Analysis of Immune Suppression in Breast Tumors and Clinical Outcomes*. This study explores the impact of immune-suppressive macrophages on breast tumors, integrating multiomics data such as RNA-seq deconvolution, cell-cell communication assays, imaging mass cytometry (IMC), and survival analysis. Follow the steps below to recreate the study’s results.

## Overview of Analysis Steps

### 1. Estimating Cell Fractions
To estimate cell fractions for bulk RNA-seq profiles, use Signature Matrices available in our previous work (https://github.com/YounessAzimzade/XML-TME-NAC-BC). We recommend using **CIBERSORTx** for cell fraction estimation. Once fractions are estimated and clinical outcome data are retrieved, you should have datasets similar to `NAC.csv` for the NAC cohort and `MBRC.csv`/`TCGA.csv` for other cohorts.

### 2. Analyzing the Role of Cell Types in Clinical Outcomes
Explore the influence of cell type frequencies on clinical outcomes using **Support Vector Machines (SSVM)** and **Random Survival Forests (RSF)**:
   - **SSVM.ipynb** and **RSF.ipynb**: Train models on cell fractions, clinical features, and outcomes to predict clinical outcomes.
   - **SHAP analysis**: Feature importance is extracted using SHAP values, yielding SHAP datasets with processed feature columns (done in SSVM.ipynb and RSF.ipynb)  
   - After SHAP analysis, **fig2.R** is used to calculate "Survival Scores."

### 3. Juxtaposing Survival and pCR Scores
Once "Survival Scores" are calculated, these scores are compared with pathological complete response (pCR) scores (from https://github.com/YounessAzimzade/XML-TME-NAC-BC) in **fig3.R**.

### 4. Classical Survival Analysis
For additional insights, run traditional survival analyses using scripts in the **Survival Analysis** folder.

### 5. Imaging Mass Cytometry (IMC) and scRNA-seq Analysis
IMC and single-cell RNA-seq data are analyzed with scripts available in their respective folders. Follow folder-specific instructions to reproduce analyses and visualizations for spatial organization and immune profiling insights.

## Repository Structure

- **data/**  
  Includes links or files for bulk RNA-seq (TCGA and MBRC), single-cell RNA-seq, and immune profiling data.

- **scripts/**  
  Contains scripts for each stage of analysis:
  - **SSVM.ipynb** and **RSF.ipynb**: Trains models to predict clinical outcomes and SHAP values.
  - **fig2.R**: Calculates Survival Scores using SHAP values.
  - **fig3.R**: Compares Survival Scores and pCR scores.
  - **IMC**: Spatial analysis of cell phenotypes with cell-cell distance calculations
    - Calculate min distance macs vs all.R
    - Cell Type Fractions.R
    - Fig 4f.R
    - Fig 4g.R
    - Fig 4g.R
    - Fig 5a.R
    - Fig 5c.R
    - Fig 5d & e.R
      
  - **scRNA-seq**:  Transfers annotation, explores cell-frequency correlations and identifies macrophage-T cell interactions via NicheNet.
    - Annotation Transfer.R
    - Fig5 g & h.R
    - Fig5 h, i, j & K.R
    - Fig5 l, m, n, o.R

## Data Access

Data sources:
- **TCGA** and **MBRC** datasets: Available on cBioPortal (https://www.cbioportal.org/).
- **Spatial omics data**:
  - Danenberg et al., 2022 (https://zenodo.org/record/5850952)
  - Wang et al., 2023 (https://zenodo.org/records/7990870)
- **scRNA-seq data**: 
  - Wu et al, 2021 (https://singlecell.broadinstitute.org/single_cell/study/SCP1039/a-single-cell-and-spatially-resolved-atlas-of-human-breast-cancers)
  - Bassez et al, 2021 (https://lambrechtslab.sites.vib.be/en/single-cell)

## Contact

For questions or support, please contact Youness Azimzade at [younessazimzade@gmail.com].
 
