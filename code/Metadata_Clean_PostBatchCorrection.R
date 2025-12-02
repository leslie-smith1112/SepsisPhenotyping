## Final Cleaning metadata; checking columns and things
# Metadata clean after batch correction -- needs the batch corrected matrix and trims metadata to samples we want
# Only disease samples

## 1.0 Read in original metadata ----
metadata <- readRDS(here::here("data","all_samples_metadata.rds"))
#metadata <- readRDS(here::here("data","metadata.rds"))


# 1.3 Make sreadRDS()# 1.3 Make sure all metadata samples are in expression matrix
dim(metadata)
colnames(metadata)
# reading in expression matrix to keep metadata that we have metadata for
expr <- readRDS(here::here("data","all_samples_expression.rds"))
expr[1:5,1:5]
dim(expr)
#trim metadata
metadata <- metadata[metadata$`Sample ID` %in% colnames(expr),]
dim(metadata)


## 2.0 Cut down metadata to only the diseases/samples we want to keep ----


# 2.1 Disease ----

meta_dat <- metadata[metadata$`Disease Simplified` %in% c("Sepsis", "Sepsis - AKI", "Sepsis - ARDs", "Sepsis - Shock"),]
table(meta_dat$`Disease Simplified`)
dim(meta_dat)

meta_dat$`Disease Simplified`[meta_dat$`Disease Simplified` == "Sepsis - ARDs"] <- "Sepsis"
meta_dat$`Disease Simplified`[meta_dat$`Disease Simplified` == "Sepsis - AKI"] <- "Sepsis"

# re-level the disease simplified column
meta_dat$`Disease Simplified` <- factor(meta_dat$`Disease Simplified`, levels = c("Sepsis", "Sepsis - Shock"))
levels(meta_dat$`Disease Simplified`)
table(meta_dat$`Disease Simplified`)
#Sepsis   Sepsis - AKI  Sepsis - ARDs Sepsis - Shock 
#2293              5             99            476


# 2.2 Get rid of Melioidosis samples ----

meta_dat <- meta_dat[!(meta_dat$`Sample Class` %in% c("Melioidosis")) ,]
table(meta_dat$`Disease Simplified`)
dim(meta_dat)
#Sepsis   Sepsis - AKI  Sepsis - ARDs Sepsis - Shock 
#2251              5             99            476 


# 2.3 Get rid of repeat samples ----

meta_dat <- meta_dat[meta_dat$`Is Repeat Sample` == "No",]
table(meta_dat$`Disease Simplified`)
dim(meta_dat)#2272   55
#Sepsis   Sepsis - AKI  Sepsis - ARDs Sepsis - Shock 
#1873              5             86            308 


# 2.4 Get rid of samples that have COVID ----

table(meta_dat$`Disease`)
meta_dat <- meta_dat[!(meta_dat$Disease %in% c("COVID Sepsis")),] # we lose 11 samples here
dim(meta_dat) #2261   55



## 3.0 Look at all values in every column and summarize values when needed (cleaning) ----


table(meta_dat$Source) #Values: GEO(777) Refine.Bio(738) SRA(553) UF(193) 
table(meta_dat$`Reanalyzed by another study`) #No(2210)  Yes (51):  
table(meta_dat[meta_dat$`Reanalyzed by another study` == "Yes",]$Dataset)
#SRP198776  expected reruns and was resolved previously, 
#GSE33118 - third party analysis we dont use any other of the datasets
#GSE33118 - is reanalyzed by another study we dont use
table(meta_dat$Tissue)
table(meta_dat$`Sample Class`)
table(meta_dat$Disease)
table(meta_dat$`Disease Simplified`)
table(meta_dat$`Disease 2`)
table(meta_dat$ARDs)


## 3.1 ARDs needs to correct capitalization and additional column to label Faheems samples ----


meta_dat$ARDs[meta_dat$ARDs == "ARDS"] <- "ARDs" # we lose 13 ARDs patients because they are repeat samples 99 -> 86

## add column for UF ARDs
meta_dat$UF_ARDs <- NA
meta_dat$UF_ARDs[meta_dat$ARDs == "ARDs"] <- "ARDs"
meta_dat$UF_ARDs[meta_dat$Dataset == "Faheem" & meta_dat$ARDs == "ARDs"] <- "UF_ARDs"
table(meta_dat$UF_ARDs)


# 3.2 fix sample missing in the sepsis/septic shock column ----
table(meta_dat$`Sepsis/Septic Shock`) #Sepsis 1897 Septic Shock 353 - we are missing one sample here

## look at sepsis and septic shock samples to see which sample is missing
septicShock <- meta_dat$`Sample ID`[meta_dat$`Sepsis/Septic Shock` == "Septic Shock"]
sepsis <- meta_dat$`Sample ID`[meta_dat$`Sepsis/Septic Shock` == "Sepsis"]
my_sampl <- meta_dat$`Sample ID`[! (meta_dat$`Sample ID` %in% c(septicShock, sepsis))] 
my_sampl 

## we see the missing sample is GSM812721 and fix it:
meta_dat$`Sepsis/Septic Shock`[meta_dat$`Sample ID` == "GSM812721"]
meta_dat$Disease[meta_dat$`Sample ID` == "GSM812721"] #this has been fixed in the metadata, but I will fix it here too since we have used this file  
meta_dat$`Sepsis/Septic Shock`[meta_dat$`Sample ID` == "GSM812721"] <- "Sepsis" # not sure what happened, was labeled as septic ards. - is marked as ARDs 
table(meta_dat$`Sepsis/Septic Shock`) # all good now

# 3.3 All these columns look good ----

table(meta_dat$Timepoint) #0 hr  0hr 12hr 16hr 24hr  2hr we have timepoints for only 1536 samples, 1079 of these are 24h
table(meta_dat$`Is Repeat Sample`)
table(meta_dat$`Low Platelets`)
table(meta_dat$Treatment)
table(meta_dat$Responder)
colnames(meta_dat)
table(meta_dat$Survival)
table(meta_dat$Dataset,meta_dat$Survival) 


# 3.4 Fix survival labels ----

## Some of the studies have survival labeled differently - describing that here:
#GSE65682 values 0 and 1, labeled as "mortalityevent28days" taking the 1 to mean death and the 0 to mean alive
#Faheem No and. Yes - labeled as "inHOspitalDeath" taking No to mean alive and yes to mean death
#GSE137340 labeled as Survived or Died
#GSE185263 labeled as Survived or Died
#GSE222393 labeled as Deceased or Survived
#GSE236713 labeled as Survived or Died
#GSE33118 negatrive or positive  labeled as survival - so taking negative as death and positive as live
#GSE95233 survivor / non-survivor

## Changing all Values to Died and Survived
#Changed to Died:   1 -> Died.    Yes -> Died.   Deceased -> Died.   negative -> Died.   Non Survivor -> Died
#Changed to Survived: 0 -> Survived.   No -> Survivied.  positive -> Survived   Survivor -> Survived
meta_dat$Survival[meta_dat$Survival == "0"] <- "Survived"
table(meta_dat$Survival)
meta_dat$Survival[meta_dat$Survival == "No"] <- "Survived"
meta_dat$Survival[meta_dat$Survival == "positive"] <- "Survived"
meta_dat$Survival[meta_dat$Survival == "Survivor"] <- "Survived"
table(meta_dat$Survival)
meta_dat$Survival[meta_dat$Survival == "1"] <- "Died"
meta_dat$Survival[meta_dat$Survival == "Deceased"] <- "Died"
meta_dat$Survival[meta_dat$Survival == "Yes"] <- "Died"
meta_dat$Survival[meta_dat$Survival == "negative"] <- "Died"
meta_dat$Survival[meta_dat$Survival == "Non Survivor"] <- "Died"
table(meta_dat$Survival)
colnames(meta_dat)


# 3.5 Make all gender columns uniform ----

table(meta_dat$Gender)
# values: F      female      Female      FEMALE           H           M        male        Male        MALE       other Transgender 
#Changing all to be Female/Male/Transgender - not sure what H and other are so leaving those as well
meta_dat$Gender[meta_dat$Gender == "F"] <- "Female"
meta_dat$Gender[meta_dat$Gender == "female"] <- "Female"
meta_dat$Gender[meta_dat$Gender == "FEMALE"] <- "Female"
meta_dat$Gender[meta_dat$Gender == "F"] <- "Female"
meta_dat$Gender[meta_dat$Gender == "M"] <- "Male"
meta_dat$Gender[meta_dat$Gender == "male"] <- "Male"
meta_dat$Gender[meta_dat$Gender == "MALE"] <- "Male"
meta_dat$Gender[meta_dat$Gender == "M"] <- "Male"
table(meta_dat$Gender)

# 3.6 Adding age summary column
meta_dat$AgeSummary <- "NoAgeReported"
meta_dat$AgeSummary[meta_dat$Age < 40] <- "Young"
meta_dat$AgeSummary[meta_dat$Age > 39 & meta_dat$Age < 65] <- "Middle-Age"
meta_dat$AgeSummary[meta_dat$Age > 64] <- "Older"



# 3.7 These all look fine: ----

table(meta_dat$`ifng il10 ratio`) # the study containing this info was eliminated for gene inclusion I believe
table(meta_dat$Ethinicity) #leaving this and race as is for now
table(meta_dat$Race)
table(meta_dat$Pathogen)
table(meta_dat$`Covid Status`) # all negatvie as it should be   
table(meta_dat$`SOFA (24hr)`)
table(meta_dat$`SOFA (24hr)`, meta_dat$Dataset) #GSE185263 and Faheem have this info
table(meta_dat$RepSOFA, meta_dat$Dataset) #only faheem
table(meta_dat$SAPSII)
table(meta_dat$SAPSII, meta_dat$Dataset) # only GSE57065
table(meta_dat$`Organ Failure Count`)
table(meta_dat$`Organ Failure Count`, meta_dat$Dataset) #only GSE63311 
table(meta_dat$`Cardiovascular Failure`, meta_dat$Dataset) #GSE63311
table(meta_dat$`Liver Failure`, meta_dat$Dataset)#GSE63311
table(meta_dat$`AKI`, meta_dat$Dataset)#GSE63311
table(meta_dat$`Coagulation Failure`, meta_dat$Dataset)#GSE63311
table(meta_dat$`Respiratory Failure`, meta_dat$Dataset)#GSE63311
table(meta_dat$`SIRS Sign`, meta_dat$Dataset)#GSE63311
table(meta_dat$`Culture Positive`, meta_dat$Dataset)#GSE63311
table(meta_dat$`Pneumonia Diagnosis`, meta_dat$Dataset) #GSE65682(cap and Hap) & GSE33118 (pneumonia)
table(meta_dat$`Thrombocytonpenia`, meta_dat$Dataset) #GSE65682

table(meta_dat$`Endotype Cohort`, meta_dat$Dataset) #GSE65682 - discovery or validation
table(meta_dat$`Endotype Class`, meta_dat$Dataset) #GSE65682(Mars) & GSE236892(hyper/hypo) & Faheem(Hypo Normo)


# 3.8 Fix Faheems endotypes because there is another study that uses hypo ----

meta_dat$`Endotype Class`[meta_dat$Dataset == "Faheem"& meta_dat$`Endotype Class` == "Normo"] <- "UF_Normo"
meta_dat$`Endotype Class`[meta_dat$Dataset == "Faheem"& meta_dat$`Endotype Class` == "Hypo"] <- "UF_Hypo"


# 3.9 These all look fine ----

table(meta_dat$`Time to Event`, meta_dat$Dataset)  #GSE65682 time to mortality
table(meta_dat$`ICU Aquired Infection`, meta_dat$Dataset)  #GSE65682
table(meta_dat$`ICU Aquired Infection Paired`, meta_dat$Dataset) #GSE65682
table(meta_dat$`Diabeties Mellitus`, meta_dat$Dataset) #GSE65682
table(meta_dat$`Abdominal Sepsis`, meta_dat$Dataset)#GSE65682
table(meta_dat$`WBC Cells per UL`, meta_dat$Dataset) #GSE66890 
table(meta_dat$`Direct Lung`, meta_dat$Dataset)#GSE66890 
table(meta_dat$`AECC ALI`, meta_dat$Dataset)#GSE66890 
table(meta_dat$`Berlin ARDs`, meta_dat$Dataset)#GSE66890 
table(meta_dat$`Shock`, meta_dat$Dataset)#GSE66890 
table(meta_dat$`Severity`, meta_dat$Dataset)#GSE69063
table(meta_dat$`UF_ARDs`, meta_dat$Dataset) 
table(meta_dat$Gender)

# final dim check -> 2261 x 56
dim(meta_dat)

## * -- save RDS -- * 
saveRDS(meta_dat,here::here("data","disease_metadata.rds")) 


