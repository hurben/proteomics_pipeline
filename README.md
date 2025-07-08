# P24-188 Proteomics Analysis

This repository contains scripts and workflows for the P24-188 proteomics data analysis project. The main components include PCA visualization, heatmap generation, linear regression modeling, volcano plot creation, and input data preparation.

## 📁 Project Structure (Example)
```
P24-188/
├── analysis/ # Snakemake pipeline (main workflow)
│ └── Snakefile # Controls the execution of all steps
├── src/ # Core analysis scripts
│ ├── PCA.r
│ ├── heatmap_v2.r
│ ├── linear_regression.r
│ ├── volcano_plot.r
│ └── prepare_input.py
```



## 🧪 Scripts Overview

### `prepare_input.py`
Prepares and formats input data (e.g., raw proteomics matrices and metadata) for downstream R scripts. Ensures column consistency and merging of metadata.

### `linear_regression.r`
Performs linear regression (or optionally linear mixed-effects modeling) to evaluate associations between protein expression and clinical/phenotypic variables. Outputs result tables for each variable.

### `PCA.r`
Generates PCA plots from the normalized protein expression matrix. Color-coded by clinical metadata such as `GENDER_CODE`, `Lookup_Race`, and binned `AGE`.

### `heatmap_v2.r`
Draws a heatmap of selected top variable/differential proteins across samples. Annotates samples with available metadata.

### `volcano_plot.r`
Creates a volcano plot based on differential expression or regression output. Allows for flexible thresholding (e.g., log2 fold-change, adjusted p-value) and feature highlighting.

## 🛠 Workflow Automation

The entire analysis is managed by a Snakemake pipeline located in the `analysis/` folder. Running the Snakefile will automatically execute all required steps in the correct order, assuming input files and environment are properly configured.

```bash

