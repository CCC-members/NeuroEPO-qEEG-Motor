#Author Fuleah A Razzaq
library(tidyverse)#load libraries
library(sjlabelled)
library(purrr)
library(readr)
library(dplyr)
library(vtable)
rm(list = ls(all.names = TRUE)) #clear memory
setwd('')#set working directory
load("")#load data
##reverse code cognition ####
cognitionTB<- cognitionTB %>% # reverse code item
  mutate(across(colour:B_A, ~ max(.x, na.rm = TRUE) - .))
##reverse code Motor
motorTB<- motorTB %>% # reverse code item
  mutate(across(speech:constancyRestTremor, ~ 4 - .))

#####set labels#####
cognitionTB<- cognitionTB  %>% 
  var_labels(
    ID = "Key column:Numeric ID extracted from raw Subject Code",
    time = "Key column: 1=baseline(T1),2= Week7(T2),3=Week30(T3),4=Week79(T4)",
    MMSE = "MMSE total Score(raw)",
    Atention = "DRS sub-score(raw)",
    Initiation = "DRS sub-score(raw)",
    Construction = "DRS sub-score(raw)",
    Conceptualization = "DRS sub-score(raw)",
    DRS_Memory = "DRS sub-score(raw)",
    DRS = "DRS total-score(raw)",              
    Arithmetic = "Working Memory WAIS III sub-score(raw)",
    digitspan = "Working Memory WAIS III sub-score(raw)",
    sequency = "Working Memory WAIS III sub-score(raw)",         
    Trial1 = "RAVLT sub-score(raw)",
    Sum = "RAVLT sub-score (Sum Trial1 to Trial5)(raw)",          
    ListB = "RAVLT sub-score(raw)",     
    recall = "RAVLT sub-score(raw)",
    Recognition = "RAVLT sub-score(raw)",
    Copy = "Rey Complex Figure sub-score(raw)",              
    Memory =  "Rey Complex Figure sub-score(raw)",
    FAB = "Frontal Assessment Battery total score(raw)",
    Phonemic = "Verbal fluency sub-score(raw)",         
    fluency = "Verbal fluency sub-score(raw)",
    Switch = "Verbal fluency sub-score(raw)",
    colour = "Stroop sub-score(reverse coded)",        
    colours_errors = "Stroop sub-score(reverse coded)",
    word = "Stroop sub-score(reverse coded)",
    word_errors = "Stroop sub-score(reverse coded)",      
    Inhibition = "Stroop sub-score(reverse coded)",
    Inhibition_errors = "Stroop sub-score(reverse coded)",
    inhibition_switch = "Stroop sub-score(reverse coded)",
    total_errors = "Stroop sub-score(reverse coded)",
    total_time = "Trail Making Test sub-score(reverse coded)",        
    B_time = "Trail Making Test sub-score(reverse coded)",           
    B_A = "Trail Making Test sub-score(reverse coded)"  
  ) 



motorTB<- motorTB  %>% 
  var_labels(
    ID = "Key column:Numeric ID extracted from raw Subject Code",
    time = "Key column: 1=baseline(T1),2= Week7(T2),3=Week30(T3),4=Week79(T4)",
    speech = "speech(0=severe, 4=normal)", 
    facialExp = "facial expression(0=severe, 4=normal)",
    rigidNeck = "rigidity neck(0=severe, 4=normal)", 
    rigidUpperL = "rigidity upper limbs left(0=severe, 4=normal)", 
    rigidUpperR = "rigidity upper limbs right(0=severe, 4=normal)",
    rigidLowerL = "rigidity lower limbs left(0=severe, 4=normal)",
    rigidLowerR = "rigidity lower limbs right(0=severe, 4=normal)",
    fingerTapL = "finger tapping left(0=severe, 4=normal)", 
    fingerTapR = "finger tapping right(0=severe, 4=normal)", 
    handMoveL = "hand movements left(0=severe, 4=normal)",
    handMoveR = "hand movements right(0=severe, 4=normal)", 
    proneSupineLH = "prone  supine left hand(0=severe, 4=normal)", 
    proneSupineRH = "prone  supine right hand(0=severe, 4=normal)",
    tapLF = "tapping left foot(0=severe, 4=normal)",
    tapRF = "tapping right foot(0=severe, 4=normal)", 
    agilityLL = "agility left leg(0=severe, 4=normal)",
    agilityRL = "agility right leg(0=severe, 4=normal)",
    chair = "Getting-up from  chair(0=severe, 4=normal)", 
    walk = "walking(0=severe, 4=normal)",
    gait = "freezing of gait(0=severe, 4=normal)",
    posturalInstability = "postural instability(0=severe, 4=normal)",
    postura = "postura(0=severe, 4=normal)", 
    movementSpontaneity = "global spontaneity of movement(0=severe, 4=normal)",
    posturalTremorLH = "postural tremor left hand(0=severe, 4=normal)", 
    posturalTremorRH = "postural tremor right hand(0=severe, 4=normal)",
    actionTremorLH = "action tremor left hand(0=severe, 4=normal)",
    actionTremorRH = "action tremor right hand(0=severe, 4=normal)", 
    restingTremorUpperL = "amplitude resting tremor left upper limbs(0=severe, 4=normal)",
    restingTremorUpperR = "amplitude resting  tremor right upper limbs(0=severe, 4=normal)", 
    restingTremorLowerL = "amplitude resting tremor left lower limb(0=severe, 4=normal)",
    restingTremorLowerR = "amplitude resting tremor right lower limb(0=severe, 4=normal)", 
    restingTremorLip = "lip/mandibular amplitude resting tremor(0=severe, 4=normal)",
    constancyRestTremor = "Constancy of rest tremor(0=severe, 4=normal)"
  )


####save data and description ####
#check directories
dir <- file.path("outReverse") 
if (!dir.exists(dir)) dir.create(dir)
dir <- file.path("outReverse","Description") 
if (!dir.exists(dir)) dir.create(dir)
dir <- file.path("outReverse","Data") 
if (!dir.exists(dir)) dir.create(dir)

vt(cognitionTB,lush=T,file="./outReverse/Description/cognition_Desc")
vt(cognitionTB,lush=T,out="csv",file="./outReverse/Description/cognitionRev_Desc.csv")
vt(subTB,lush=T,file="./outReverse/Description/SubjectAttr_Desc")
vt(subTB,lush=T,out="csv",file="./outReverse/Description/SubjectAttr_Desc.csv")
vt(subTimeTB,lush=T,file="./outReverse/Description/SubjectTimeAttr_Desc")
vt(subTimeTB,lush=T,out="csv",file="./outReverse/Description/SubjectTimeAttr_Desc.csv")
vt(motorTB,lush=T,file="./outReverse/Description/MotorAttr_Desc")
vt(motorTB,lush=T,out="csv",file="./outReverse/Description/MotorRevAttr_Desc.csv")

save(cognitionTB,subTB,subTimeTB,motorTB,file="./outReverse/Data/tidyEPOdataRev.Rdata")
write.table(cognitionTB , file = "./outReverse/Data/cognitionTBRev.csv", sep=",", row.names=FALSE)
write.table(motorTB , file = "./outReverse/Data/motorTBRev.csv", sep=",", row.names=FALSE)
write.table(subTB , file = "./outReverse/Data/SubjectTB.csv", sep=",", row.names=FALSE)
write.table(subTimeTB , file = "./outReverse/Data/SubjectTimeTB.csv", sep=",", row.names=FALSE)

