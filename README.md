# NeuroEPO and Parkinson’s Disease: Causal Analysis of Motor Performance

## Overview  
This repository contains MATLAB and R scripts used for the secondary analysis of a **randomized controlled trial (RCT)** investigating the causal effect of **NeuroEPO** on **motor performance** in **Parkinson’s Disease (PD)**. The analysis employs **marginal structural models (MSMs)** and **mediation analysis** to explore causal relationships between treatment, motor function, and brain activity as reflected in **quantitative EEG (qEEG) measures**.

## Contents  

### 1. Scripts for Data Processing & Analysis  
- **MATLAB scripts:**  
  - `XYZ_to_AAL_BrainNet_coord.m` – Converts XYZ coordinates to **AAL BrainNet Viewer format** for visualization.  
  - `nii_gen_4_BrainNet.m` – Generates `.nii` files for **BrainNet Viewer**, mapping **CCA scores** onto the cortex.  

- **R scripts:**  
  - `mediation_analysis_Commented.R` – Performs **mediation analysis**, testing whether qEEG mediates the effect of **NeuroEPO** on motor performance.  
  - `msmIPWTcommented.R` – Implements **marginal structural models (MSMs)** to estimate the **causal effect** of NeuroEPO.  

## Key Findings  
- NeuroEPO significantly improves **motor performance**, particularly **bradykinesia and rigidity** symptoms.  
- The **EEG signal mediates 93% of the observed motor improvements**, suggesting a **direct neural mechanism** for NeuroEPO.  
- The study employs **Factor Analysis** to extract a **latent motor variable** summarizing motor function based on **UPDRS-III items**.  
- **Marginal Structural Models (MSMs)** estimate causal effects, while **mediation models** confirm EEG as a mediator.  

## Variable Name Equivalence Table  
This table maps the variable names used in the **article** and **code**, along with a short description of each:

| **Article Variable Name**     | **Code Variable Name**     | **Description** |
|------------------------------|--------------------------|----------------|
| Motor latent variable        | `Lamda_motor`   | Factor-analyzed motor score derived from UPDRS-III items. |
| Treatment Group              | `Doseind`                | Binary variable indicating treatment with **NeuroEPO (1) or Placebo (0)**. |
| Motor assessments (UPDRS-III) | `speech`, `rigidNeck`, `fingerTapL`, etc. | Individual motor test items from UPDRS-III. |
| EEG Canonical Component (qEEG mediator) | `KqEEG` | First **Canonical Correlation Analysis (CCA) latent factor** from qEEG. |
| Motor Canonical Component    | `Kmotor`     | First **CCA latent factor** from motor scores. |
| Causal effect estimate       | `msm`                   | **Marginal Structural Model (MSM)** estimating the causal effect of NeuroEPO. |
| Inverse probability weights (IPTW) | `siptw`           | Adjusted weights used in **MSM causal modeling**. |
| Side of onset of symptoms    | `side`                   | Indicates whether Parkinson’s symptoms started on the left or right side. |
| Disease severity             | `severity`               | Parkinson’s Disease severity based on **Hoehn & Yahr scale**. |
| Confounders                  | `progression`, `age`, `handness` | Additional covariates adjusted for in the MSM analysis. |

## Data & Dependencies  
### Data Sources  
- The study used **UPDRS-III motor assessment data** and **quantitative EEG (qEEG) recordings**.  
- Original datasets are available in **OpenNeuro and OSF repositories**:  
  - [NeuroEPO EEG data (OpenNeuro)](https://openneuro.org/datasets/ds003194/versions/1.0.0)  
  - [Placebo EEG data (OpenNeuro)](https://openneuro.org/datasets/ds003195/versions/1.0.0)  
  - [Tidy format datasets (OSF)](https://doi.org/10.17605/OSF.IO/M8SJP)  

### Required Software & Libraries  
#### MATLAB (for .nii file generation & coordinate transformation)  
- **Toolboxes Required:**  
  - `Tools for NIfTI and Analyze image` (for `.nii` file handling)  
  - `BrainNet Viewer` (for visualization)  

#### R (for causal & mediation analysis)  
- **Required Packages:**  
  - `tidyverse`, `dplyr`, `readr` – Data manipulation  
  - `mediation` – Mediation analysis  
  - `lmerTest`, `lme4` – Mixed effects models  
  - `ipw` – Marginal Structural Models  
  - `psych` – Factor analysis  
  - `EValue` – Sensitivity analysis  

## Usage Instructions  
1. **Preprocess Data:**  
   - Ensure **.Rdata files** (e.g., `tidyEPOdata.Rdata`) are available in the working directory.  

2. **Run Causal Models:**  
   - Execute `msmIPWTcommented.R` to **estimate the causal effect** of NeuroEPO.  
   - Execute `mediation_analysis_Commented.R` to analyze **EEG mediation effects**.  

3. **Visualize Results:**  
   - Run MATLAB scripts to **generate .nii files** for BrainNet Viewer visualization.  

4. **Interpret Findings:**  
   - Review outputs and load `.nii` files into **BrainNet Viewer** for 3D brain visualizations.  

## Authors & Contributors  
This project was developed by:  
- **Maria L. Bringas-Vega, PhD**  
- **Fuleah A. Razzaq, PhD**  
- **Liu Shengnan, MSc**  
- **Carlos Lopez Naranjo, PhD**  
- **Pedro A. Valdes-Sosa, PhD**  

For inquiries, contact:  
**pedro.valdes@neuroinformatics-collaboratory.org**  

## Citing This Work  
If you use this code or data, please cite:  
> Bringas-Vega, M. L., Razzaq, F. A., Shengnan, L., et al. "Causal analysis of the effect of NeuroEPO on motor performance in Parkinson’s Disease: Secondary analysis of a randomized controlled study."  

## License  
This repository is released under the **MIT License**. See `LICENSE` for details.  
