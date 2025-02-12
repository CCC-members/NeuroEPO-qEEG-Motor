# Author: Fuleah A. Razzaq
# Modified by Carlos Lopez
library(tidyverse) # Load necessary libraries
library(readr)
library(dplyr)

rm(list = ls(all.names = TRUE)) # Clear memory
load("tidyEPOdata.Rdata") # Load dataset

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

Y=as.matrix(subset(data,select=(speech:dyskinesia))) # UPDRS-III motor outcome variables (Methods: "Motor assessment")

library(psych)
my_fa <- fa(r = Y, nfactors = 3)
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

library(geepack)
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
