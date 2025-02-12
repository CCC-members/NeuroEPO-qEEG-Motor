# NeuroEPO and Parkinson’s Disease: Causal Analysis of Motor Performance
Overview
This repository contains MATLAB and R scripts used for the secondary analysis of a randomized controlled trial (RCT) investigating the causal effect of NeuroEPO on motor performance in Parkinson’s Disease (PD). The analysis employs marginal structural models (MSMs) and mediation analysis to explore causal relationships between treatment, motor function, and brain activity as reflected in quantitative EEG (qEEG) measures.
Contents
This repository includes:
1. Scripts for Data Processing & Analysis
•	MATLAB scripts:
o	XYZ_to_AAL_BrainNet_coord.m – Converts XYZ coordinates to AAL BrainNet Viewer format for visualization.
o	nii_gen_4_BrainNet.m – Generates .nii files for BrainNet Viewer, mapping CCA scores onto the cortex.
•	R scripts:
o	mediation_analysis_Commented.R – Performs mediation analysis, testing whether qEEG mediates the effect of NeuroEPO on motor performance.
o	msmIPWTcommented.R – Implements marginal structural models (MSMs) to estimate the causal effect of NeuroEPO.
o	ReverseCodeMotor.R – Prepares and reverse-codes cognition and motor data for factor analysis.
Key Findings
•	NeuroEPO significantly improves motor performance, particularly bradykinesia and rigidity symptoms.
•	The EEG signal mediates 93% of the observed motor improvements, suggesting a direct neural mechanism for NeuroEPO.
•	The study employs Factor Analysis to extract a latent motor variable summarizing motor function based on UPDRS-III items.
•	Marginal Structural Models (MSMs) estimate causal effects, while mediation models confirm EEG as a mediator.
Data & Dependencies
Data Sources
•	The study used UPDRS-III motor assessment data and quantitative EEG (qEEG) recordings.
•	Original datasets are available in OpenNeuro and OSF repositories:
o	NeuroEPO EEG data (OpenNeuro)
o	Placebo EEG data (OpenNeuro)
o	Tidy format datasets (OSF)
Required Software & Libraries
MATLAB (for .nii file generation & coordinate transformation)
•	Toolboxes Required:
o	Tools for NIfTI and Analyze image (for .nii file handling)
o	BrainNet Viewer (for visualization)
R (for causal & mediation analysis)
•	Required Packages:
o	tidyverse, dplyr, readr – Data manipulation
o	mediation – Mediation analysis
o	lmerTest, lme4 – Mixed effects models
o	ipw – Marginal Structural Models
o	psych – Factor analysis
o	EValue – Sensitivity analysis
Usage Instructions
1.	Preprocess Data:
o	Use ReverseCodeMotor.R to prepare motor.
o	Ensure .Rdata files (e.g., tidyEPOdata.Rdata) are available in the working directory.
2.	Run Causal Models:
o	Execute msmIPWTcommented.R to estimate the causal effect of NeuroEPO.
o	Execute mediation_analysis_Commented.R to analyze EEG mediation effects.
3.	Visualize Results:
o	Run MATLAB scripts to generate .nii files for BrainNet Viewer visualization.
4.	Interpret Findings:
o	Review outputs and load .nii files into BrainNet Viewer for 3D brain visualizations.
Authors & Contributors
This project was developed by:
•	Maria L. Bringas-Vega, PhD
•	Fuleah A. Razzaq, PhD
•	Liu Shengnan, MSc
•	Carlos Lopez Naranjo, MSc
•	Pedro A. Valdes-Sosa, PhD
For inquiries, contact:
pedro.valdes@neuroinformatics-collaboratory.org
Citing This Work
If you use this code or data, please cite:
Bringas-Vega, M. L., Razzaq, F. A., Shengnan, L., et al. "Causal analysis of the effect of NeuroEPO on motor performance in Parkinson’s Disease: Secondary analysis of a randomized controlled study."

