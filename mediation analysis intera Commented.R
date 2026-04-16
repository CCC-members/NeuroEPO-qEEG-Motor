# Author: Fuleah A. Razzaq
# Modified by Carlos Lopez

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

data <- motorTB %>%
  inner_join(eegTB, by = c("ID", "time")) %>%
  # In the tracked EEG tables, the post-intervention EEG visit is coded as time = 3.
  filter(time %in% c(1, 3)) %>%
  na.omit()

Y <- as.matrix(data %>% dplyr::select(all_of(selected_motor_items)))
X <- as.matrix(data %>% dplyr::select(EEG.X1:EEG.X158956))

cca.out <- whitening::scca(X, Y, scale = TRUE)
wx1 <- extract_component_loadings(cca.out$WX, component = 1L, expected_length = ncol(X), label = "X")
wy1 <- extract_component_loadings(cca.out$WY, component = 1L, expected_length = ncol(Y), label = "Y")

resCCA <- data.frame(
  ID = data$ID,
  time = data$time,
  KqEEG = as.vector(scale(X, center = TRUE, scale = FALSE) %*% wx1),
  Kmotor = as.vector(scale(Y, center = TRUE, scale = FALSE) %*% wy1)
)

resCCA <- subTB %>% inner_join(resCCA, by = "ID")
resCCA <- subTimeTB %>% inner_join(resCCA, by = c("ID", "time"))
resCCA <- resCCA %>%
  filter(complete.cases(ID, time, KqEEG, Kmotor, Dose, progression, handness, side, severity, age)) %>%
  mutate(
    ID = factor(ID),
    side = factor(side),
    handness = factor(handness)
  )

model.m.lmer <- as.formula("KqEEG ~ 1 + Dose + progression + handness + side + severity + age + (1 | ID)")
model.o.lmer <- as.formula("Kmotor ~ 1 + KqEEG * Dose + progression + handness + side + severity + age + (1 | ID)")

fit.m <- lme4::lmer(model.m.lmer, data = resCCA, REML = FALSE)
fit.o <- lme4::lmer(model.o.lmer, data = resCCA, REML = FALSE)

results <- mediation::mediate(
  fit.m,
  fit.o,
  treat = "Dose",
  mediator = "KqEEG",
  sims = 5000,
  na.action = "na.omit"
)

summary(results)
print(summary(fit.m))
print(summary(fit.o))

single_cca_score2 <- wx1 * cca.out$lambda[1]
writeMat(file.path(root_dir, "single_cca_scores.mat"), single_cca_score2 = single_cca_score2)

saveRDS(
  list(
    selected_motor_items = selected_motor_items,
    wx1 = wx1,
    wy1 = wy1,
    fit.m = fit.m,
    fit.o = fit.o,
    mediation = results,
    resCCA = resCCA
  ),
  file.path(root_dir, "mediation_interaction_resCCA.rds")
)

print("single_cca_scores.mat has been successfully created!")