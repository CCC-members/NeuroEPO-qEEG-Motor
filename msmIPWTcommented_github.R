# GitHub-targeted copy of msmIPWTcommented.R
# Keeps the original script untouched while making the workflow repo-portable
# and the latent motor scoring explicit for manuscript-aligned upload.

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

script_path <- get_script_path()
root_dir <- dirname(script_path)

.libPaths(c(file.path(root_dir, "Rlibs"), .libPaths()))

library(tidyverse) # Load necessary libraries
library(readr)
library(dplyr)

# Load datasets (Methods: "Data Acquisition")
# load(file.path(root_dir, "outCloud/newData/tidyEPOdataOn.Rdata")) # Load EPO dataset
load(file.path(root_dir, "tidyEPOdata.Rdata")) # Load dataset


# Rename variables while keeping original read names intact
colnames(subTB)[7] = "side" # Side of symptom onset (Methods: "Confounder variables")

subgr=subTimeTB
subgr$group=NA
subgr$group[subgr$time==2 & subgr$Dose==0]="Placebo" # Treatment variable (Methods: "Participants")
subgr$group[subgr$time==2 & subgr$Dose==5]="NeuroEPO"
subgr<-na.omit(subgr)%>%dplyr::select(ID,group)
subTB<-subgr %>%  
  inner_join(subTB, by = ('ID'))

# Merge motor and cognitive test data
data<-cognitionTB %>%  
  inner_join(motorTB, by = c('ID','time'))
data=data[data$time %in% c(1,2,3),]

motor_item_names <- names(motorTB)[match("speech", names(motorTB)):match("dyskinesia", names(motorTB))]
Y=as.matrix(data[, motor_item_names]) # 33 OFF-state MDS-UPDRS Part III motor items (Methods: "Motor assessment")

library(psych)
my_fa <- psych::fa(
  Y,
  nfactors = 3,
  fm = "minres",
  rotate = "oblimin",
  scores = "regression"
)

if (is.null(my_fa$scores) || nrow(my_fa$scores) != nrow(data)) {
  stop("psych::fa did not return factor scores for all longitudinal rows.")
}

resY=my_fa$scores # Factor analysis to extract latent motor variable (Methods: "Exploring evidence for a causal effect")

####### Combine data #####
resCCA=data.frame(ID=data$ID,time=data$time,resY)

#### Covariates ####
resCCA<-subTB %>%  
  inner_join(resCCA, by = ('ID'))
resCCA<-subTimeTB %>%  
  inner_join(resCCA, by = c('ID','time'))

resCCA$Doseind=0
resCCA$Doseind[resCCA$Dose>0]=1
resCCA$Doseind=as.integer(resCCA$Doseind) # Binary treatment indicator (Methods: "Statistical Analysis")
resCCA$severity=as.integer(resCCA$severity) # Disease severity (Methods: "Confounder variables")
resCCA$handeness=as.factor(resCCA$handness) # Handedness (Methods: "Confounder variables")
resCCA$side=as.factor(resCCA$side) # Side of onset of symptoms (Methods: "Confounder variables")
resCCA$days=resCCA$age-resCCA$initage # Derived variable: age at study - age at onset

# Define previous treatment and motor assessments (Methods: "Marginal Structural Models")
resCCA$prevX=0
resCCA$prevY=resCCA$prevY2=resCCA$prevY3=NA
for (id in 1:26 )
{
  resCCA$prevX[resCCA$ID==id & resCCA$time %in% c(3,4)]=resCCA$Doseind[resCCA$ID==id & resCCA$time==2]
  resCCA$prevY[resCCA$ID==id & resCCA$time %in% c(1,2,3,4)]=resCCA$MR1[resCCA$ID==id & resCCA$time==1]# MR1 represents the motor latent variable derived from UPDRS-III (Methods: "Exploring evidence for a causal effect of NeuroEPO on the motor outcome")
  resCCA$prevY2[resCCA$ID==id & resCCA$time %in% c(1,2,3,4)]=resCCA$MR2[resCCA$ID==id & resCCA$time==1]
  resCCA$prevY3[resCCA$ID==id & resCCA$time %in% c(1,2,3,4)]=resCCA$MR3[resCCA$ID==id & resCCA$time==1]
}

dataCausal=data.frame(resCCA)
dataCausal <- dataCausal %>% rename(Lamda_motor = MR1) # Lamda_motor is MR1

library(ipw)
library(EValue)
# Compute inverse probability weights for MSM (Methods: "Marginal Structural Models")
siptw <- ipwtm(exposure = Doseind,timevar=time,
               family = "binomial",  link="logit",
               numerator= ~1 + prevX,
               denominator =  ~  prevX +prevY+progression+age+severity+handeness+side, # Adjusted confounders (Methods: "Confounder variables")
               id = ID,  type = "all",  data = dataCausal, corstr="AR1")

quantile(siptw$ipw.weights) # Check distribution of weights (Methods: "Marginal Structural Models")

dataCausal$siptw <- siptw$ipw.weights
# MSM estimation using linear mixed models (Methods: "Statistical Analysis")
msm=nlme::lme(Lamda_motor~Doseind, data=dataCausal, 
              weights =~siptw,
              random=~1|ID)
summary(msm)
s=summary(msm)

# Sensitivity analysis for unmeasured confounders (Methods: "Sensitivity Analysis")
beta_lme = s$tTable["Doseind", "Value"]
se_lme = s$tTable["Doseind", "Std.Error"]
sd_lme = s$sigma

e_lme = evalues.OLS(est = beta_lme, 
                    se = se_lme, 
                    sd = sd_lme, 
                    delta = 1)

library(boot) # Load bootstrapping package
library(nlme)  # Load the package that contains fixef()


# Function to estimate MSM for bootstrapped samples
msm_bootstrap <- function(data, indices) {
  data_boot <- data[indices, ]  # Resample data
  
  # Compute IPTW weights
  siptw <- ipwtm(exposure = Doseind, timevar = time,
                 family = "binomial", link = "logit",
                 numerator = ~1 + prevX,
                 denominator = ~prevX + prevY + progression + age + severity + handness + side,
                 id = ID, type = "all", data = data_boot, corstr="AR1")
  
  data_boot$siptw <- siptw$ipw.weights
  
  # MSM estimation using linear mixed-effects model
  msm_model <- nlme::lme(Lamda_motor ~ Doseind, data = data_boot, 
                         weights = ~siptw, random = ~1 | ID)
  
  return(fixef(msm_model)["Doseind"]) # Return ATE estimate
}

# Run bootstrap with 1000 iterations
set.seed(123)  # For reproducibility
boot_results <- boot(dataCausal, msm_bootstrap, R = 1000)

# Compute 95% confidence intervals
boot_ci <- boot.ci(boot_results, type = c("perc", "bca"))

# Print results
print(boot_results)
print(boot_ci)