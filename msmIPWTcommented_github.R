# =============================================================================
# msmIPWTcommented_github.R
# =============================================================================
# Purpose:
#   Estimate the causal effect of NeuroEPO (vs placebo) on the latent motor
#   outcome (Lambda_motor) using a Marginal Structural Model (MSM) fitted via
#   Inverse Probability of Treatment Weighting (IPTW). The latent motor
#   variable is extracted from 33 OFF-state MDS-UPDRS Part III items using
#   three-factor oblimin-rotated factor analysis.
#   A bootstrap confidence interval is computed for the average treatment
#   effect (ATE) and an E-value is reported for sensitivity to unmeasured
#   confounding (Methods: "Marginal Structural Models", "Sensitivity Analysis").
#
# Inputs:
#   tidyEPOdata.Rdata -- longitudinal clinical/demographic EPO trial data
#     Exposes: subTB, subTimeTB, motorTB, cognitionTB
#
# Outputs (console):
#   MSM summary, E-value, bootstrap CI
#
# Author: Fuleah A. Razzaq
# Modified by Carlos Lopez
# =============================================================================

# ---------------------------------------------------------------------------
# Helper: resolve the path of the current script so data files are located
# relative to the script, making the workflow repo-portable.
# ---------------------------------------------------------------------------
get_script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  matches <- grep(file_arg, args, value = TRUE)

  if (length(matches) == 0) {
    stop("Could not determine script path from commandArgs().")
  }

  normalizePath(sub(file_arg, "", matches[1]), winslash = "/", mustWork = TRUE)
}

script_path <- get_script_path()
root_dir <- dirname(script_path)

# Point R to the project-local library so no system-wide installation is needed.
.libPaths(c(file.path(root_dir, "Rlibs"), .libPaths()))

# ---------------------------------------------------------------------------
# Load packages
# ---------------------------------------------------------------------------
library(tidyverse) # data wrangling and plotting utilities
library(readr)     # CSV/flat-file I/O
library(dplyr)     # tidy data manipulation

# ---------------------------------------------------------------------------
# Load dataset (Methods: "Data Acquisition")
# tidyEPOdata.Rdata exposes: subTB (subject-level), subTimeTB (visit-level),
# motorTB (MDS-UPDRS Part III items per visit), cognitionTB (cognitive items).
# ---------------------------------------------------------------------------
# load(file.path(root_dir, "outCloud/newData/tidyEPOdataOn.Rdata")) # Load EPO dataset
load(file.path(root_dir, "tidyEPOdata.Rdata")) # Load dataset

# Ensure 'side' (side of motor symptom onset) is correctly named in subTB
# (Methods: "Confounder variables").
colnames(subTB)[7] = "side" # Side of symptom onset (Methods: "Confounder variables")

# ---------------------------------------------------------------------------
# Derive treatment group assignment from the second visit (time = 2).
# Dose == 0 -> Placebo; Dose == 5 -> NeuroEPO
# (Methods: "Participants")
# ---------------------------------------------------------------------------
subgr=subTimeTB
subgr$group=NA
subgr$group[subgr$time==2 & subgr$Dose==0]="Placebo" # Treatment variable (Methods: "Participants")
subgr$group[subgr$time==2 & subgr$Dose==5]="NeuroEPO"
subgr<-na.omit(subgr)%>%dplyr::select(ID,group)
subTB<-subgr %>%  
  inner_join(subTB, by = ('ID'))

# ---------------------------------------------------------------------------
# Merge motor and cognitive items into one longitudinal data frame, keeping
# visits 1, 2, and 3 (baseline, mid-point, and 6-month follow-up in the MSM)
# (Methods: "Motor assessment", "Exploring evidence for a causal effect")
# ---------------------------------------------------------------------------
data<-cognitionTB %>%  
  inner_join(motorTB, by = c('ID','time'))
data=data[data$time %in% c(1,2,3),]

# ---------------------------------------------------------------------------
# Latent motor factor extraction
# (Methods: "Exploring evidence for a causal effect of NeuroEPO on the motor outcome")
#
# All 33 OFF-state MDS-UPDRS Part III motor items are assembled into matrix Y.
# Three latent factors are extracted with oblimin rotation; the first factor
# (MR1) is taken as the latent motor variable Lambda_motor.
# Scoring is done by regression on the observed item scores.
# ---------------------------------------------------------------------------
motor_item_names <- names(motorTB)[match("speech", names(motorTB)):match("dyskinesia", names(motorTB))]
Y=as.matrix(data[, motor_item_names]) # 33 OFF-state MDS-UPDRS Part III motor items (Methods: "Motor assessment")

library(psych)
my_fa <- psych::fa(
  Y,
  nfactors = 3,     # Three factors chosen to capture tremor, rigidity, and axial motor dimensions
  fm = "minres",    # Minimum residual factoring (robust to non-normality of ordinal items)
  rotate = "oblimin", # Oblique rotation allowing correlated factors
  scores = "regression" # Factor scores estimated by regression for best linear prediction
)

# Abort if factor scores were not returned for every longitudinal row.
if (is.null(my_fa$scores) || nrow(my_fa$scores) != nrow(data)) {
  stop("psych::fa did not return factor scores for all longitudinal rows.")
}

resY=my_fa$scores # Factor analysis to extract latent motor variable (Methods: "Exploring evidence for a causal effect")

# ---------------------------------------------------------------------------
# Assemble analysis data frame: combine factor scores with covariates.
# ---------------------------------------------------------------------------
####### Combine data #####
resCCA=data.frame(ID=data$ID,time=data$time,resY) # resY columns: MR1 (Lambda_motor), MR2, MR3

#### Covariates ####
resCCA<-subTB %>%  
  inner_join(resCCA, by = ('ID'))   # Append subject-level variables (group, Dose)
resCCA<-subTimeTB %>%  
  inner_join(resCCA, by = c('ID','time')) # Append visit-level variables (age, severity, progression)

# Encode treatment as a binary integer (0 = placebo, 1 = NeuroEPO)
# required by ipwtm() (Methods: "Marginal Structural Models").
resCCA$Doseind=0
resCCA$Doseind[resCCA$Dose>0]=1
resCCA$Doseind=as.integer(resCCA$Doseind) # Binary treatment indicator (Methods: "Statistical Analysis")
resCCA$severity=as.integer(resCCA$severity)   # Disease severity at baseline (Methods: "Confounder variables")
resCCA$handeness=as.factor(resCCA$handness)   # Handedness as factor (Methods: "Confounder variables")
resCCA$side=as.factor(resCCA$side)             # Side of symptom onset as factor (Methods: "Confounder variables")
resCCA$days=resCCA$age-resCCA$initage         # Derived: years from symptom onset to study entry

# ---------------------------------------------------------------------------
# Build lagged treatment and outcome history variables required by the MSM.
# For each participant, prevX records the treatment received at the previous
# visit, and prevY / prevY2 / prevY3 record the baseline latent factor scores.
# These are used in both the IPTW numerator/denominator models and as
# time-varying confounders in the MSM (Methods: "Marginal Structural Models").
# ---------------------------------------------------------------------------
# Define previous treatment and motor assessments (Methods: "Marginal Structural Models")
resCCA$prevX=0
resCCA$prevY=resCCA$prevY2=resCCA$prevY3=NA
for (id in 1:26 )
{
  # prevX: treatment indicator at the preceding visit (used in IPTW denominator)
  resCCA$prevX[resCCA$ID==id & resCCA$time %in% c(3,4)]=resCCA$Doseind[resCCA$ID==id & resCCA$time==2]
  # prevY / prevY2 / prevY3: baseline factor scores used as time-stable confounders
  resCCA$prevY[resCCA$ID==id & resCCA$time %in% c(1,2,3,4)]=resCCA$MR1[resCCA$ID==id & resCCA$time==1]# MR1 represents the motor latent variable derived from UPDRS-III (Methods: "Exploring evidence for a causal effect of NeuroEPO on the motor outcome")
  resCCA$prevY2[resCCA$ID==id & resCCA$time %in% c(1,2,3,4)]=resCCA$MR2[resCCA$ID==id & resCCA$time==1]
  resCCA$prevY3[resCCA$ID==id & resCCA$time %in% c(1,2,3,4)]=resCCA$MR3[resCCA$ID==id & resCCA$time==1]
}

dataCausal=data.frame(resCCA)
dataCausal <- dataCausal %>% rename(Lamda_motor = MR1) # Lamda_motor is MR1 (first latent motor factor)

library(ipw)    # ipwtm(): inverse probability of treatment weighting
library(EValue) # evalues.OLS(): E-value for sensitivity analysis

# ---------------------------------------------------------------------------
# IPTW computation
# Stabilised IPTW weights are computed using a logistic model for treatment
# assignment. The denominator model adjusts for all baseline and time-varying
# confounders; the numerator model adjusts for treatment history only
# (Methods: "Marginal Structural Models").
# ---------------------------------------------------------------------------
# Compute inverse probability weights for MSM (Methods: "Marginal Structural Models")
siptw <- ipwtm(exposure = Doseind,timevar=time,
               family = "binomial",  link="logit",
               numerator= ~1 + prevX,
               denominator =  ~  prevX +prevY+progression+age+severity+handeness+side, # Adjusted confounders (Methods: "Confounder variables")
               id = ID,  type = "all",  data = dataCausal, corstr="AR1")

quantile(siptw$ipw.weights) # Check distribution of weights -- large values may indicate positivity violations (Methods: "Marginal Structural Models")

dataCausal$siptw <- siptw$ipw.weights

# ---------------------------------------------------------------------------
# MSM estimation
# The weighted linear mixed-effects model yields the average treatment effect
# (ATE) of NeuroEPO on Lambda_motor under the MSM, marginalised over the
# distribution of time-varying confounders (Methods: "Statistical Analysis").
# ---------------------------------------------------------------------------
# MSM estimation using linear mixed models (Methods: "Statistical Analysis")
msm=nlme::lme(Lamda_motor~Doseind, data=dataCausal, 
              weights =~siptw,    # Incorporate stabilised IPTW weights
              random=~1|ID)       # Random intercept per participant
summary(msm)
s=summary(msm)

# ---------------------------------------------------------------------------
# E-value sensitivity analysis
# Quantifies the minimum strength of association an unmeasured confounder
# would need to have with both treatment and outcome to fully explain the
# observed ATE (Methods: "Sensitivity Analysis").
# ---------------------------------------------------------------------------
# Sensitivity analysis for unmeasured confounders (Methods: "Sensitivity Analysis")
beta_lme = s$tTable["Doseind", "Value"]    # Estimated ATE
se_lme = s$tTable["Doseind", "Std.Error"]  # Standard error of ATE
sd_lme = s$sigma                            # Residual SD of the MSM

e_lme = evalues.OLS(est = beta_lme, 
                    se = se_lme, 
                    sd = sd_lme, 
                    delta = 1) # delta = 1: one-unit change in the unmeasured confounder

library(boot)  # Load bootstrapping package -- boot(), boot.ci()
library(nlme)  # Load nlme for lme() and fixef()

# ---------------------------------------------------------------------------
# Bootstrap confidence interval for the ATE
# Re-estimates IPTW weights and the MSM within each bootstrap replicate to
# propagate uncertainty from the weight model into the final CI
# (Methods: "Statistical Analysis").
# ---------------------------------------------------------------------------

# Function to estimate MSM for bootstrapped samples
msm_bootstrap <- function(data, indices) {
  data_boot <- data[indices, ]  # Resample data (rows sampled with replacement)
  
  # Compute IPTW weights on the resampled data
  siptw <- ipwtm(exposure = Doseind, timevar = time,
                 family = "binomial", link = "logit",
                 numerator = ~1 + prevX,
                 denominator = ~prevX + prevY + progression + age + severity + handness + side,
                 id = ID, type = "all", data = data_boot, corstr="AR1")
  
  data_boot$siptw <- siptw$ipw.weights
  
  # Fit weighted MSM on the bootstrap sample
  msm_model <- nlme::lme(Lamda_motor ~ Doseind, data = data_boot, 
                         weights = ~siptw, random = ~1 | ID)
  
  return(fixef(msm_model)["Doseind"]) # Return bootstrap ATE estimate
}

# Run bootstrap with 1000 iterations for a stable percentile and BCa CI
set.seed(123)  # For reproducibility
boot_results <- boot(dataCausal, msm_bootstrap, R = 1000)

# Compute 95% confidence intervals using percentile and bias-corrected-accelerated methods
boot_ci <- boot.ci(boot_results, type = c("perc", "bca"))

# Print results
print(boot_results)
print(boot_ci)