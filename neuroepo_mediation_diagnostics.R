############################################################
# FULL DIAGNOSTIC + ROBUSTNESS SCRIPT FOR THE EEG-MOTOR
# MEDIATION ANALYSIS REPORTED IN THE MANUSCRIPT.
#
# This script rebuilds resCCA directly from tidyEPOdata.Rdata and
# tidyEEGVARETA.Rdata, using the screened OFF-state motor items that
# enter the sparse CCA and the original nomenclature KqEEG / Kmotor.
############################################################

get_script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  matches <- grep(file_arg, args, value = TRUE)

  if (length(matches) == 0) {
    stop("Could not determine script path from commandArgs().")
  }

  normalizePath(sub(file_arg, "", matches[1]), winslash = "/", mustWork = TRUE)
}

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

.libPaths(c(file.path(root_dir, "Rlibs"), .libPaths()))

req_pkgs <- c(
  "dplyr",
  "whitening",
  "lme4",
  "lmerTest",
  "DHARMa",
  "performance",
  "clubSandwich",
  "influence.ME",
  "mediation",
  "splines",
  "lmtest"
)

to_install <- req_pkgs[!req_pkgs %in% installed.packages()[, "Package"]]
if (length(to_install) > 0) install.packages(to_install, repos = "https://cran.rstudio.com/")

suppressPackageStartupMessages({
  library(dplyr)
  library(whitening)
  library(lme4)
  library(lmerTest)
  library(DHARMa)
  library(performance)
  library(clubSandwich)
  library(influence.ME)
  library(mediation)
  library(splines)
  library(lmtest)
})

load(file.path(root_dir, "tidyEPOdata.Rdata"))
load(file.path(root_dir, "tidyEEGVARETA.Rdata"))

if (!"side" %in% names(subTB) && ncol(subTB) >= 7) {
  colnames(subTB)[7] <- "side"
}

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

missing_motor_items <- setdiff(selected_motor_items, names(motorTB))
if (length(missing_motor_items) > 0) {
  stop(sprintf("Missing screened motor items in motorTB: %s", paste(missing_motor_items, collapse = ", ")))
}

analysis_data <- motorTB %>%
  inner_join(eegTB, by = c("ID", "time")) %>%
  # In the tracked EEG tables, the post-intervention EEG visit is coded as time = 3.
  filter(time %in% c(1, 3)) %>%
  na.omit()

X <- as.matrix(analysis_data %>% dplyr::select(EEG.X1:EEG.X158956))
Y <- as.matrix(analysis_data %>% dplyr::select(all_of(selected_motor_items)))
X <- as.matrix(analysis_data %>% dplyr::select(EEG.X1:EEG.X158956))

cca.out <- whitening::scca(X, Y, scale = TRUE)
wx1 <- extract_component_loadings(cca.out$WX, component = 1L, expected_length = ncol(X), label = "X")
wy1 <- extract_component_loadings(cca.out$WY, component = 1L, expected_length = ncol(Y), label = "Y")

resCCA <- data.frame(
  ID = analysis_data$ID,
  time = analysis_data$time,
  KqEEG = as.vector(scale(X, center = TRUE, scale = FALSE) %*% wx1),
  Kmotor = as.vector(scale(Y, center = TRUE, scale = FALSE) %*% wy1)
)

resCCA <- subTB %>%
  inner_join(resCCA, by = "ID")
resCCA <- subTimeTB %>%
  inner_join(resCCA, by = c("ID", "time")) %>%
  filter(complete.cases(ID, time, KqEEG, Kmotor, Dose, progression, handness, side, severity, age)) %>%
  mutate(
    ID = factor(ID),
    side = factor(side),
    handness = factor(handness)
  )

############################
# 1) CHECK THAT VARIABLES EXIST
############################

print(names(resCCA))
vars_needed <- c("KqEEG", "Kmotor", "Dose", "progression", "handness", "side", "severity", "age", "ID")
print(vars_needed %in% names(resCCA))

############################
# 2) FIT MAIN MEDIATOR AND OUTCOME MODELS
############################

fit.m <- lme4::lmer(
  KqEEG ~ 1 + Dose + progression + handness + side + severity + age + (1 | ID),
  data = resCCA,
  REML = FALSE
)

fit.o <- lme4::lmer(
  Kmotor ~ 1 + KqEEG + Dose + progression + handness + side + severity + age + (1 | ID),
  data = resCCA,
  REML = FALSE
)

############################
# 3) SUMMARIES OF BOTH MODELS
############################

cat("\n===== MEDIATOR MODEL SUMMARY =====\n")
print(summary(fit.m))

cat("\n===== OUTCOME MODEL SUMMARY =====\n")
print(summary(fit.o))

############################
# 4) BASIC DIAGNOSTICS: MEDIATOR MODEL
############################

if (interactive()) {
  plot(fit.m)
  qqnorm(residuals(fit.m), main = "Mediator model residuals")
  qqline(residuals(fit.m))
}
print(shapiro.test(residuals(fit.m)))

############################
# 5) BASIC DIAGNOSTICS: OUTCOME MODEL
############################

if (interactive()) {
  plot(fit.o)
  qqnorm(residuals(fit.o), main = "Outcome model residuals")
  qqline(residuals(fit.o))
}
print(shapiro.test(residuals(fit.o)))

############################
# 6) DHARMa DIAGNOSTICS FOR MIXED MODELS
############################

set.seed(123)
sim_m <- simulateResiduals(fittedModel = fit.m, n = 1000)
if (interactive()) plot(sim_m)
print(testDispersion(sim_m))
print(testOutliers(sim_m))

sim_o <- simulateResiduals(fittedModel = fit.o, n = 1000)
if (interactive()) plot(sim_o)
print(testDispersion(sim_o))
print(testOutliers(sim_o))

############################
# 7) COLLINEARITY CHECK IN OUTCOME MODEL
############################

cat("\n===== COLLINEARITY CHECK =====\n")
print(check_collinearity(fit.o))

############################
# 8) NONLINEARITY CHECK: QUADRATIC TERM
############################

fit.o.quad <- lme4::lmer(
  Kmotor ~ 1 + KqEEG + I(KqEEG^2) + Dose + progression + handness + side + severity + age + (1 | ID),
  data = resCCA,
  REML = FALSE
)

cat("\n===== QUADRATIC NONLINEARITY CHECK =====\n")
print(summary(fit.o.quad))
print(anova(fit.o, fit.o.quad))

############################
# 9) NONLINEARITY CHECK: SPLINE MODEL
############################

fit.o.spline <- lme4::lmer(
  Kmotor ~ 1 + ns(KqEEG, df = 3) + Dose + progression + handness + side + severity + age + (1 | ID),
  data = resCCA,
  REML = FALSE
)

cat("\n===== SPLINE NONLINEARITY CHECK =====\n")
print(summary(fit.o.spline))
print(anova(fit.o, fit.o.spline))

############################
# 10) TREATMENT-MEDIATOR INTERACTION CHECK
############################

fit.o.int <- lme4::lmer(
  Kmotor ~ 1 + KqEEG * Dose + progression + handness + side + severity + age + (1 | ID),
  data = resCCA,
  REML = FALSE
)

cat("\n===== INTERACTION CHECK =====\n")
print(summary(fit.o.int))
print(anova(fit.o, fit.o.int))

############################
# 11) CLUSTER-ROBUST SEs FOR OUTCOME MODEL
############################

vcov_cr2_o <- vcovCR(fit.o, cluster = resCCA$ID, type = "CR2")
robust_o <- coef_test(fit.o, vcov = vcov_cr2_o, test = "Satterthwaite")

cat("\n===== ROBUST SEs: OUTCOME MODEL =====\n")
print(robust_o)

############################
# 12) CLUSTER-ROBUST SEs FOR INTERACTION MODEL
############################

vcov_cr2_o_int <- vcovCR(fit.o.int, cluster = resCCA$ID, type = "CR2")
robust_o_int <- coef_test(fit.o.int, vcov = vcov_cr2_o_int, test = "Satterthwaite")

cat("\n===== ROBUST SEs: INTERACTION MODEL =====\n")
print(robust_o_int)

############################
# 13) INFLUENCE DIAGNOSTICS BY ID
############################

infl_m <- influence(fit.m, group = "ID")
cat("\n===== INFLUENCE: MEDIATOR MODEL =====\n")
print(summary(infl_m))
if (interactive()) plot(infl_m, which = "cook")

infl_o <- influence(fit.o, group = "ID")
cat("\n===== INFLUENCE: OUTCOME MODEL =====\n")
print(summary(infl_o))
if (interactive()) plot(infl_o, which = "cook")

############################
# 14) QUASI-BAYESIAN MEDIATION
############################

set.seed(123)
med.qb <- mediate(
  model.m = fit.m,
  model.y = fit.o,
  treat = "Dose",
  mediator = "KqEEG",
  sims = 5000
)

cat("\n===== MEDIATION: QUASI-BAYESIAN =====\n")
print(summary(med.qb))

############################
# 15) BOOTSTRAP MEDIATION
############################

cat("\n===== MEDIATION: BOOTSTRAP =====\n")
med.boot <- NULL
cat(
  "Skipped: mediation::mediate does not support boot = TRUE for lme4::lmer models. ",
  "The quasi-Bayesian mediation fit above matches the main mediation workflow.\n",
  sep = ""
)

############################
# 16) OPTIONAL: MEDIATION WITH INTERACTION MODEL
############################

set.seed(123)
med.int <- mediate(
  model.m = fit.m,
  model.y = fit.o.int,
  treat = "Dose",
  mediator = "KqEEG",
  sims = 5000
)

cat("\n===== MEDIATION WITH INTERACTION MODEL =====\n")
print(summary(med.int))

############################
# 17) SAVE RESULTS
############################

all_checks <- list(
  manuscript_context = list(
    eeg_time_points = c(1L, 3L),
    selected_motor_items = selected_motor_items,
    final_reported_model = "interaction-adjusted mediation model"
  ),
  resCCA = resCCA,
  cca = list(fit = cca.out, wx1 = wx1, wy1 = wy1),
  fit.m = fit.m,
  fit.o = fit.o,
  fit.o.quad = fit.o.quad,
  fit.o.spline = fit.o.spline,
  fit.o.int = fit.o.int,
  sim_m = sim_m,
  sim_o = sim_o,
  robust_o = robust_o,
  robust_o_int = robust_o_int,
  infl_m = infl_m,
  infl_o = infl_o,
  med.qb = med.qb,
  med.boot = med.boot,
  med.int = med.int
)

saveRDS(all_checks, file.path(root_dir, "mediation_diagnostics_resCCA.rds"))