# NeuroEPO and Parkinson’s Disease: Causal Analysis of Motor Performance

## Overview

This repository contains MATLAB and R scripts for the secondary analysis of a **randomized controlled trial (RCT)** investigating the causal effect of **NeuroEPO** on **motor performance** in **Parkinson's Disease (PD)**. The analysis employs **Marginal Structural Models (MSMs)** with **Inverse Probability of Treatment Weighting (IPTW)**, **sparse Canonical Correlation Analysis (sCCA)**, and **causal mediation analysis** to characterize the relationship between treatment, latent motor function (λ_motor), and brain activity captured by **quantitative EEG (qEEG)** source spectra.

**Study:** *Causal analysis of the effect of NeuroEPO on motor performance in Parkinson's Disease: Secondary analysis of a randomized controlled study*
**Running title:** *Causal effect of NeuroEPO in Parkinson's Disease*

---

## Contents

### 1. Scripts for Data Processing & Analysis

#### R scripts

| Script | Purpose |
|---|---|
| `msmIPWTcommented_github.R` | **MSM/IPTW causal analysis.** Extracts latent motor variable (λ_motor) via 3-factor oblimin FA on 33 OFF-state MDS-UPDRS Part III items; computes stabilized IPTW weights; fits weighted linear mixed-effects MSM; reports ATE, E-value, and 1000-replicate bootstrap BCa CI. |
| `mediation analysis intera Commented.R` | **Causal mediation analysis (interaction-adjusted).** Runs sCCA to derive κ_qEEG and κ_motor; fits mediator model (κ_qEEG ~ group_neuroepo + covariates) and interaction-adjusted outcome model (κ_motor ~ κ_qEEG × group_neuroepo + covariates); estimates ACME, ADE, and proportion mediated via 5000 Monte Carlo simulations. |
| `ExportLocfdrCcaWeights.R` | **CCA weight thresholding.** Applies local false discovery rate (locfdr) to the raw CCA X-block weight vector; exports significant weights, FDR values, and mask to `.mat` for NIfTI generation. |
| `reviewer_delta_lambda_vs_updrs_analysis.R` | **Reviewer validation.** Computes participant-level Δλ_motor (T4−T1) vs. ΔUPDRS-III OFF total (T4−T1); reports Pearson r with 95% CI and Spearman ρ with bootstrap 95% CI; generates scatter plots with regression CI band, colored by treatment group. |
| `neuroepo_mediation_diagnostics_github.R` | **Mediation diagnostics.** Self-contained; rebuilds the sCCA pipeline from `.Rdata` files; runs influence analysis, robust SEs, and DHARMa residual checks. |
| `scca_first_pair_inference_nuisance_adjusted.R` | **sCCA inference.** Permutation-based inference for the first CCA pair after nuisance-variable adjustment. |

#### MATLAB scripts

| Script | Purpose |
|---|---|
| `CreateNiisLocfdrFigureMaxMin.m` | Reads `cca_wx_locfdr_component0_1.mat`; reshapes 3244-voxel × 49-frequency weights by frequency band (Δ, Θ, α, β); writes signed, max, and min NIfTI maps to `CCA_locfdr_figure_maxmin_maps0_1/`. |
| `PlotBrainNetLocfdrDisplayBatch4ViewsBipolar.m` | Reads NIfTI maps; renders 4-view BrainNet figures (left sagittal, right sagittal, axial top, axial bottom) with a symmetric blue-white-red-yellow bipolar colormap. |
| `J2img3d_nors_sg3244.m` | NIfTI writer utility for the 3244-voxel VARETA source space. |
| `CreateNiisPosNeg.m` | Generates positive/negative NIfTI maps for BrainNet visualization. |

---

## Key Findings

### Causal Effect of NeuroEPO on Motor Function (MSM/IPTW)

The MSM estimated an **average treatment effect (ATE) of −0.49** on the latent motor variable λ_motor (p = 0.014; 95% CI: −0.66 to −0.21). Because λ_motor is aligned so that higher values indicate worse motor function, the negative sign indicates **motor improvement with NeuroEPO**. Stabilized IPTW weights ranged 0.5–3.7.

Bootstrap validation (1000 replicates, BCa method): **95% CI −0.91 to −0.12** (bias = 0.003, SE = 0.20).

**E-value = 3.47** — an unmeasured confounder would require an association of at least 3.47 with both treatment and outcome to explain away the observed effect, indicating robustness to unmeasured confounding.

> *Note:* Conventional UPDRS-III total scores showed within-group improvements of approximately −9.7 to −9.9 points in both arms, with a non-significant between-group difference (ΔΔ = +0.18, p = 0.964), highlighting that the latent-variable causal model captures signal not visible in raw group comparisons.

### EEG as a Causal Mediator (Interaction-Adjusted Mediation)

The interaction-adjusted mediation model (κ_qEEG × randomized group term in outcome model) yielded:

| Parameter | Estimate | p-value |
|---|---|---|
| **ACME** (qEEG-mediated effect) | **−0.107** | **0.0036** |
| **ADE** (direct treatment effect) | **−0.029** | **0.0012** |
| **Proportion mediated** | **~79%** | — |

The EEG latent variable (κ_qEEG) mediates approximately **79% of the total causal effect** of NeuroEPO on motor function, consistent with a direct neural mechanism for the drug's action.

### CCA Motor-EEG Coupling

The top MDS-UPDRS Part III items loaded on the first sCCA motor canonical variate (κ_motor) are dominated by **axial and appendicular bradykinesia/rigidity** (CCA loadings 0.62–0.74): Rigid Neck (0.74), Facial Expression (0.73), Finger Tap R (0.70), Prone-Supine R Hand (0.68), Rigid Lower R (0.66), Speech (0.64), Hand Move R (0.62).

### Study Sample

- 26 participants randomized (NeuroEPO n = 15, Placebo n = 11).
- 1 placebo participant excluded for EEG artifacts → **final analytical n = 25** (NeuroEPO = 15, Placebo = 10).
- Baseline MDS-UPDRS III OFF-state total median: 44 [IQR 37–50]; disease duration median: 5 years; LED median: 788 mg.

---

## Variable Name Equivalence Table

| Article Name | Code Name | Description |
|---|---|---|
| Latent motor variable (λ_motor) | `Lamda_motor` / `MR1` | First factor from 3-factor oblimin FA on 33 OFF-state MDS-UPDRS Part III items. Higher = worse motor function. |
| EEG latent variable (κ_qEEG) | `KqEEG` | First sCCA canonical variate from qEEG source spectra; the mediator. |
| Motor canonical variate (κ_motor) | `Kmotor` | First sCCA canonical variate from screened motor items; the outcome in mediation. |
| Group | `group`; model indicator: `group_neuroepo` | Randomized study group. `group` is coded as Placebo or NeuroEPO; `group_neuroepo` is coded 0 = Placebo and 1 = NeuroEPO for causal and mediation models. |
| Average Treatment Effect | `beta_lme` (MSM) | Estimated by weighted LME; negative = motor improvement. |
| ACME | `results$d.avg` | Average Causal Mediation Effect from `mediation::mediate()`. |
| ADE | `results$z.avg` | Average Direct Effect from `mediation::mediate()`. |
| Stabilized IPTW weights | `siptw` | From `ipwtm()` with AR1 correlation structure. |
| E-value | `e_lme` | From `evalues.OLS()`; quantifies unmeasured-confounding robustness. |
| MDS-UPDRS Part III items | `speech`, `facialExp`, `rigidNeck`, `fingerTapL/R`, etc. | 33 OFF-state motor items used for λ_motor; 18 screened items used for sCCA Y-block. |
| Disease severity (HY) | `severity` | Hoehn & Yahr stage; integer covariate. |
| Side of onset | `side` | Left/right; factor covariate. |
| Handedness | `handness` / `handeness` | Factor covariate. |
| Additional confounders | `progression`, `age` | Time-varying covariates in IPTW denominator model. |
| EEG source spectra | `EEG.X1:EEG.X158956` | 3244 voxels × 49 frequency bins (0.78–19.1 Hz, Δf = 0.39 Hz) from VARETA inverse solution. |
| locfdr significant mask | `wx_sig_mask` | Binary; 1 where local FDR < threshold for CCA weight. |

The exposure variable was the randomized group, coded as placebo group = 0 and NeuroEPO group = 1 for causal and mediation models.

---

## Data & Dependencies

### Data Sources

- **Clinical + motor data:** `tidyEPOdata.Rdata` (exposes `subTB`, `subTimeTB`, `motorTB`, `cognitionTB`)
- **EEG source spectra:** `tidyEEGVARETA.Rdata` (exposes `eegTB`; 3244-voxel × 49-frequency VARETA source spectra)
- **Pre-computed CCA outputs:** `cca_wx_locfdr_component0_1.mat`, `single_cca_scores.mat`
- Original raw EEG datasets on OpenNeuro:
  - [NeuroEPO EEG (ds003194)](https://openneuro.org/datasets/ds003194/versions/1.0.0)
  - [Placebo EEG (ds003195)](https://openneuro.org/datasets/ds003195/versions/1.0.0)

### Required R Packages

| Package | Purpose |
|---|---|
| `tidyverse`, `dplyr`, `readr` | Data wrangling |
| `psych` | Factor analysis (λ_motor) |
| `whitening` | sCCA via `whitening::scca()` |
| `lme4`, `lmerTest` | Mixed-effects models (mediation) |
| `nlme` | Mixed-effects models (MSM) |
| `ipw` | Stabilized IPTW via `ipwtm()` |
| `mediation` | Causal mediation (ACME/ADE) |
| `EValue` | E-value sensitivity analysis |
| `boot` | Bootstrap BCa confidence intervals |
| `locfdr` | Local FDR thresholding of CCA weights |
| `R.matlab` | Read/write `.mat` files |
| `ggplot2` | Visualization |

### Required MATLAB Toolboxes / Software

- **BrainNet Viewer** — brain surface visualization
- **Tools for NIfTI and Analyze image** — `.nii` file I/O
- **VARETA** — Variable Resolution Electromagnetic Tomography (EEG source spectra)
- **EEGLAB** — EEG preprocessing and artifact rejection

---

## Usage Instructions

### Step 1 — MSM/IPTW Causal Analysis
```r
# Estimate causal effect of NeuroEPO on λ_motor
Rscript msmIPWTcommented_github.R
```
Outputs ATE (−0.49, p = 0.014), bootstrap BCa CI, and E-value (3.47) to console.

### Step 2 — Causal Mediation Analysis (with Group × Mediator Interaction)
```r
# Estimate ACME, ADE, and proportion mediated by κ_qEEG
Rscript "mediation analysis intera Commented.R"
```
Outputs ACME (−0.107, p = 0.0036), ADE (−0.029, p = 0.0012), proportion mediated (~79%), and saves `mediation_interaction_resCCA.rds`.

### Step 3 — Export locfdr-Thresholded CCA Weights
```r
# Requires cca.out from whitening::scca(); run inside mediation script or separately
source("ExportLocfdrCcaWeights.R")
# Outputs: cca_wx_locfdr_component0_1.mat
```

### Step 4 — Generate NIfTI Brain Maps (MATLAB)
```matlab
% Run after cca_wx_locfdr_component0_1.mat exists
CreateNiisLocfdrFigureMaxMin
% Outputs: CCA_locfdr_figure_maxmin_maps0_1/*.nii
```

### Step 5 — Render 4-View BrainNet Figures (MATLAB)
```matlab
% Run after NIfTI maps are generated
PlotBrainNetLocfdrDisplayBatch4ViewsBipolar
% Outputs: BrainNet_Locfdr_Display_4Views_Output0_1_BipolarOnly/**/*.png
```

### Step 6 — Reviewer Delta Analysis (optional)
```r
Rscript reviewer_delta_lambda_vs_updrs_analysis.R
# Outputs: scatter plots and Pearson/Spearman CIs in reviewer_delta_lambda_vs_updrs_output/
```

---

## Authors & Contributors

| Author | Affiliation |
|---|---|
| **Maria L. Bringas-Vega, PhD** *(shared first)* | UESTC, Chengdu; CIREN, Havana |
| **Carlos Lopez Naranjo, PhD** *(shared first)* | UESTC, Chengdu |
| Fuleah A. Razzaq, PhD | UESTC, Chengdu |
| Liu Shengnan, MSc | UESTC, Chengdu |
| Marite Garcia Llano, PhD | CIREN, Havana |
| Daniel Amaro Gonzalez, PhD | CIM, Havana |
| Lilia Morales Chacon, PhD | UESTC; CIREN; Universidad de La Rioja |
| Jorge Bosch Bayard, PhD | Montreal Neurological Institute; UESTC |
| Trinidad Virues Alba, BSc | Cuban Neuroscience Center |
| Usama Riaz, PhD | UESTC, Chengdu |
| Sambhu Pandi, MSc | UESTC, Chengdu |
| Marjan Jahanshahi, PhD | UCL Queen Square Institute of Neurology; UESTC |
| **Ivonne Pedroso Ibañez, MD** *(shared corresponding)* | CIREN, Havana |
| **Pedro A. Valdes-Sosa, PhD** *(shared corresponding)* | UESTC, Chengdu; Cuban Neuroscience Center |

For inquiries: **pedro.valdes@neuroinformatics-collaboratory.org**

---

## Citing This Work

If you use this code or data, please cite:
> Bringas-Vega, M. L., Lopez Naranjo, C., Razzaq, F. A., et al. "Causal analysis of the effect of NeuroEPO on motor performance in Parkinson's Disease: Secondary analysis of a randomized controlled study."

---

## License

This repository is released under the **MIT License**. See `LICENSE` for details.
