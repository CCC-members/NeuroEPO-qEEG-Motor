# =============================================================================
# mediation analysis intera Commented.R
# =============================================================================
# Purpose:
#   Estimate the causal mediation effect of qEEG (KqEEG) on the latent motor
#   outcome (Kmotor) through NeuroEPO randomized group, incorporating a
#   NeuroEPO-group-by-mediator interaction term (KqEEG * group_neuroepo) in the outcome model.
#   This extends the base mediation model to allow the effect of KqEEG on
#   Kmotor to vary by randomized group, in line with reviewer suggestions and
#   manuscript reporting (Methods: "Statistical Analysis").
#
# Inputs:
#   tidyEPOdata.Rdata  -- longitudinal clinical/demographic EPO trial data
#   tidyEEGVARETA.Rdata -- preprocessed qEEG feature matrix (EEG.X1:EEG.X158956)
#
# Outputs:
#   mediation_interaction_resCCA.rds -- serialised mediation object + CCA data
#   single_cca_scores.mat            -- first CCA X-canonical score (MATLAB format)
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

# ---------------------------------------------------------------------------
# Helper: safely extract the weight vector for one CCA component from
# whitening::scca output, handling both row-major and column-major storage.
# ---------------------------------------------------------------------------
extract_component_loadings <- function(weights, component, expected_length, label) {
  if (is.null(dim(weights))) {
    if (length(weights) != expected_length) {
      stop(sprintf("Unexpected %s loading length: expected %d, got %d", label, expected_length, length(weights)))
    }
    return(as.numeric(weights))
  }

  if (ncol(weights) == expected_length) {
    return(as.numeric(weights[component, , drop = TRUE]))
  }

  if (nrow(weights) == expected_length) {
    return(as.numeric(weights[, component, drop = TRUE]))
  }

  stop(sprintf(
    "Could not determine %s loading orientation from dimensions %s",
    label,
    paste(dim(weights), collapse = " x ")
  ))
}

script_path <- get_script_path()
root_dir <- dirname(script_path)

# Point R to the project-local library so no system-wide installation is needed.
.libPaths(c(file.path(root_dir, "Rlibs"), .libPaths()))

# ---------------------------------------------------------------------------
# Package management: install any missing packages automatically.
# ---------------------------------------------------------------------------
req_pkgs <- c("tidyverse", "readr", "dplyr", "whitening", "mediation", "lmerTest", "lme4", "R.matlab")
to_install <- req_pkgs[!req_pkgs %in% installed.packages()[, "Package"]]
if (length(to_install) > 0) install.packages(to_install, repos = "https://cran.rstudio.com/")

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(dplyr)
  library(whitening)
  library(mediation)
  library(lmerTest)
  library(lme4)
  library(R.matlab)
})

# ---------------------------------------------------------------------------
# Load data (Methods: "Data Acquisition")
#   tidyEPOdata.Rdata   -- contains subTB, subTimeTB, motorTB, cognitionTB
#   tidyEEGVARETA.Rdata -- contains eegTB with one row per participant-visit
# ---------------------------------------------------------------------------
load(file.path(root_dir, "tidyEPOdata.Rdata"))
load(file.path(root_dir, "tidyEEGVARETA.Rdata"))

# Ensure 'side' (side of symptom onset) is correctly named in the subject table
# (Methods: "Confounder variables").
if (!"side" %in% names(subTB) && ncol(subTB) >= 7) {
  colnames(subTB)[7] <- "side"
}

# Derive randomized group assignment from the active-treatment visit.
# The legacy source field is converted once per participant at load time.
group_lookup <- subTimeTB %>%
  filter(time == 2, !is.na(Dose)) %>%
  transmute(
    ID,
    group = factor(
      if_else(Dose > 0, "NeuroEPO", "Placebo"),
      levels = c("Placebo", "NeuroEPO")
    ),
    group_neuroepo = as.integer(group == "NeuroEPO")
  ) %>%
  distinct(ID, .keep_all = TRUE)

# ---------------------------------------------------------------------------
# Motor items used as the Y block in sparse CCA.
# These 18 MDS-UPDRS Part III items were selected via factor-analytic
# screening to represent the motor construct assessed in the OFF state
# (Methods: "Motor assessment").
# ---------------------------------------------------------------------------
selected_motor_items <- c(
  "speech",
  "facialExp",
  "rigidNeck",
  "rigidUpperL",
  "rigidUpperR",
  "fingerTapL",
  "fingerTapR",
  "handMoveL",
  "handMoveR",
  "proneSupineLH",
  "proneSupineRH",
  "tapLF",
  "tapRF",
  "chair",
  "posturalTremorLH",
  "posturalTremorRH",
  "actionTremorLH",
  "actionTremorRH"
)

# Sanity check: abort early if any expected column is absent.
missing_motor_items <- setdiff(selected_motor_items, names(motorTB))
if (length(missing_motor_items) > 0) {
  stop(sprintf("Missing screened motor items in motorTB: %s", paste(missing_motor_items, collapse = ", ")))
}

# ---------------------------------------------------------------------------
# Merge motor and EEG tables, retaining only the two EEG-recorded visits:
#   time = 1 (manuscript T1, baseline) and time = 3 (manuscript T4, 6-month).
# Note: The EEG data are available only at T1 and T4; within the tracked
# eegTB the post-intervention visit is stored as time = 3, not time = 4.
# (Methods: "EEG recording and processing")
# ---------------------------------------------------------------------------
data <- motorTB %>%
  inner_join(eegTB, by = c("ID", "time")) %>%
  # In the tracked EEG tables, the post-intervention EEG visit is coded as time = 3.
  filter(time %in% c(1, 3)) %>%
  na.omit()

# Y block: OFF-state motor items (18 selected MDS-UPDRS Part III items).
Y <- as.matrix(data %>% dplyr::select(all_of(selected_motor_items)))
# X block: high-dimensional qEEG feature matrix (158,956 spectral features).
X <- as.matrix(data %>% dplyr::select(EEG.X1:EEG.X158956))

# ---------------------------------------------------------------------------
# Sparse canonical correlation analysis (sCCA) to identify the maximally
# correlated linear combination of qEEG features (KqEEG) and motor items
# (Kmotor), as described in Methods: "Sparse Canonical Correlation Analysis".
# whitening::scca applies a whitening pre-transform before CCA to improve
# numerical stability with high-dimensional X.
# ---------------------------------------------------------------------------
cca.out <- whitening::scca(X, Y, scale = TRUE)

# Extract canonical weight vectors for the first (most correlated) component.
wx1 <- extract_component_loadings(cca.out$WX, component = 1L, expected_length = ncol(X), label = "X")
wy1 <- extract_component_loadings(cca.out$WY, component = 1L, expected_length = ncol(Y), label = "Y")

# Project the mean-centred data onto the canonical directions to obtain
# participant-visit canonical scores:
#   KqEEG  -- latent qEEG composite (mediator)
#   Kmotor -- latent motor composite (outcome)
resCCA <- data.frame(
  ID = data$ID,
  time = data$time,
  KqEEG = as.vector(scale(X, center = TRUE, scale = FALSE) %*% wx1),
  Kmotor = as.vector(scale(Y, center = TRUE, scale = FALSE) %*% wy1)
)

# Append subject-level covariates, randomized group, and visit-level covariates
# (progression, handness, side, severity, age) from the EPO data tables
# (Methods: "Confounder variables").
resCCA <- subTB %>% inner_join(resCCA, by = "ID")
resCCA <- group_lookup %>% inner_join(resCCA, by = "ID")
resCCA <- subTimeTB %>% inner_join(resCCA, by = c("ID", "time"))

# Retain only complete cases on all analysis variables and cast grouping
# variables to factors for correct treatment in the mixed-effects models.
resCCA <- resCCA %>%
  filter(complete.cases(ID, time, KqEEG, Kmotor, group_neuroepo, progression, handness, side, severity, age)) %>%
  mutate(
    ID = factor(ID),      # Random effect grouping variable (participant ID)
    side = factor(side),  # Side of motor symptom onset
    handness = factor(handness), # Handedness
    group = factor(group, levels = c("Placebo", "NeuroEPO")),
    group_neuroepo = as.integer(group_neuroepo)
  )

# ---------------------------------------------------------------------------
# Causal mediation model with randomized-group-by-mediator interaction
# (Methods: "Statistical Analysis", "Mediation Analysis")
#
# Mediator model (M-model):
#   KqEEG ~ group_neuroepo + covariates + (1 | ID)
#   -- Tests whether NeuroEPO shifts the qEEG mediator relative to placebo.
#
# Outcome model (O-model) -- interaction-adjusted:
#   Kmotor ~ KqEEG * group_neuroepo + covariates + (1 | ID)
#   -- The KqEEG:group_neuroepo interaction allows the qEEG-to-motor path to differ
#      between NeuroEPO and placebo arms, relaxing the no-interaction
#      assumption of the standard Baron-Kenny framework.
#   -- Under this model, mediation::mediate() returns
#      ACME (average causal mediation effect) and ADE (average direct effect)
#      integrated over the observed randomized group distribution.
# ---------------------------------------------------------------------------
model.m.lmer <- as.formula("KqEEG ~ 1 + group_neuroepo + progression + handness + side + severity + age + (1 | ID)")
model.o.lmer <- as.formula("Kmotor ~ 1 + KqEEG * group_neuroepo + progression + handness + side + severity + age + (1 | ID)")

# Fit both mixed-effects models using ML (REML = FALSE) for likelihood-based
# comparison and compatibility with mediation::mediate().
fit.m <- lme4::lmer(model.m.lmer, data = resCCA, REML = FALSE)
fit.o <- lme4::lmer(model.o.lmer, data = resCCA, REML = FALSE)

# Estimate ACME, ADE, and total effect via quasi-Bayesian Monte Carlo
# simulation (sims = 5000) (Methods: "Mediation Analysis").
results <- mediation::mediate(
  fit.m,
  fit.o,
  treat = "group_neuroepo",
  mediator = "KqEEG",
  sims = 5000,
  na.action = "na.omit"
)

summary(results)        # ACME, ADE, total effect and proportion mediated
print(summary(fit.m))  # Mediator model coefficient table
print(summary(fit.o))  # Outcome model coefficient table (including KqEEG:group_neuroepo)

# ---------------------------------------------------------------------------
# Export results
# ---------------------------------------------------------------------------

# Scale the first canonical X-weight by the corresponding canonical correlation
# (lambda[1]) and save to MATLAB format for downstream visualisation of the
# qEEG topographic pattern (Methods: "EEG source imaging").
single_cca_score2 <- wx1 * cca.out$lambda[1]
writeMat(file.path(root_dir, "single_cca_scores.mat"), single_cca_score2 = single_cca_score2)

# Serialise the full analysis object (CCA weights, model fits, mediation output,
# and the analysis data frame) for audit and downstream sensitivity checks.
saveRDS(
  list(
    selected_motor_items = selected_motor_items, # 18 screened MDS-UPDRS items used in Y block
    wx1 = wx1,          # First CCA X-canonical weight vector (qEEG)
    wy1 = wy1,          # First CCA Y-canonical weight vector (motor)
    fit.m = fit.m,      # Fitted mediator model
    fit.o = fit.o,      # Fitted outcome model (with KqEEG:group_neuroepo interaction)
    mediation = results, # mediation::mediate() output (ACME, ADE, total effect)
    resCCA = resCCA     # Analysis data frame (canonical scores + covariates)
  ),
  file.path(root_dir, "mediation_interaction_resCCA.rds")
)

print("single_cca_scores.mat has been successfully created!")
