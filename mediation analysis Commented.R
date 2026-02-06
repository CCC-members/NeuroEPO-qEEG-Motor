# Author: Fuleah A. Razzaq
# Modified by Carlos Lopez
.libPaths("d:/OneDrive/Downloads/NeuroEpoGitRepo/NeuroEPO-qEEG-Motor/Rlibs")
# Load necessary libraries
if (!require("tidyverse")) install.packages("tidyverse", repos = "https://cran.rstudio.com/")
if (!require("readr")) install.packages("readr", repos = "https://cran.rstudio.com/")
if (!require("dplyr")) install.packages("dplyr", repos = "https://cran.rstudio.com/")
if (!require("whitening")) install.packages("whitening", repos = "https://cran.rstudio.com/")
if (!require("mediation")) install.packages("mediation", repos = "https://cran.rstudio.com/")
if (!require("lmerTest")) install.packages("lmerTest", repos = "https://cran.rstudio.com/")
if (!require("lme4")) install.packages("lme4", repos = "https://cran.rstudio.com/")

library(tidyverse) # Load necessary libraries
library(readr)
library(dplyr)
library(whitening)
library(mediation)
library(lmerTest)

rm(list = ls(all.names = TRUE)) # Clear memory

# Load datasets (Methods: "Data Acquisition")
load("tidyEPOdata.Rdata") # Load EPO dataset
load("tidyEEGVARETA.Rdata") # Load EEG dataset

# Merge cognition, motor, and EEG data (Methods: "Data Processing")
data <- motorTB %>%  
  inner_join(eegTB, by = c('ID', 'time'))
data <- na.omit(data) # Remove missing values

# Rename side of onset variable to match terminology (Methods: "Confounder Variables")
colnames(subTB)[7] = "side"

# Define matrices for mediation analysis (Methods: "Statistical Analysis")
Y = as.matrix(subset(data, select = (speech:actionTremorRH))) # Outcome variables
X = as.matrix(subset(data, select = (EEG.X1:EEG.X158956))) # EEG predictor variables

# Canonical correlation analysis (Methods: "Canonical Correlation Analysis")
cca.out = scca(X, Y, scale = TRUE)

# Compute canonical components (Methods: "Feature Extraction")
CCAX = tcrossprod(scale(X, center = TRUE, scale = FALSE), cca.out$WX)
CCAY = tcrossprod(scale(Y, center = TRUE, scale = FALSE), cca.out$WY)

# Create dataset for mediation analysis
resCCA = data.frame(ID = data$ID, time = data$time, CCAX = CCAX, CCAY = CCAY)

# Merge covariates into dataset
resCCA <- subTB %>% inner_join(resCCA, by = ('ID'))
resCCA <- subTimeTB %>% inner_join(resCCA, by = c('ID', 'time'))
resCCA <- resCCA %>% rename(KqEEG = CCAX.1) # CCAX.1 is KqEEG
resCCA <- resCCA %>% rename(Kmotor = CCAY.1) #CCAY.1 is Kmotor

# Define mediation models (Methods: "Mediation Analysis")
model.m.lmer = as.formula(paste0("KqEEG ~ 1 + Dose + progression + handness + side + severity + age + (1 | ID)"))  
model.o.lmer = as.formula(paste0("Kmotor ~ 1 + KqEEG + Dose + progression + handness + side + severity + age + (1 | ID)")) 

# # Fit linear mixed-effects models (Methods: "Statistical Modeling")
fit.m = lme4::lmer(model.m.lmer, data = resCCA)
fit.o = lme4::lmer(model.o.lmer, data = resCCA)

# Just for doing the ANOVA, not used with the mediation, since it does not support lmerTest::lmer
# fit.m = lmerTest::lmer(model.m.lmer, data = resCCA)
# fit.o = lmerTest::lmer(model.o.lmer, data = resCCA)

# Perform mediation analysis (Methods: "Mediation Analysis")
results = mediation::mediate(fit.m, fit.o, treat = 'Dose', mediator = "KqEEG", sims = 1000, na.action = "na.omit")
summary(results)

# Compute ANOVA for model evaluation
# anova(fit.m)  # ANOVA for the mediator model
# anova(fit.o)  # ANOVA for the outcome model

# Extract p-values from models
# coef(summary(fit.m))[ , "Pr(>|t|)"]  # p-values for mediator model
# coef(summary(fit.o))[ , "Pr(>|t|)"]  # p-values for outcome model


# Ensure R.matlab package is installed
if (!require("R.matlab")) install.packages("R.matlab", dependencies = TRUE)
library(R.matlab)


# Compute the single CCA score vector (First Canonical Component)
single_cca_score2 = cca.out$WX[1, ] * cca.out$lambda[1]

# Save as a .mat file
writeMat("single_cca_scores.mat", single_cca_score2 = single_cca_score2)

# Confirm the file is saved
print("single_cca_scores.mat has been successfully created!")
