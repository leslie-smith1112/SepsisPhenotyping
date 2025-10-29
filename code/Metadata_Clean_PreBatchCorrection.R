# Metadata cleaning pre-batch correction ----
# For cleaning up some of the disease stuff in the metadata to look at during batch correction

## 1.0 Read in original metadata ----
metadata <- readr::read_csv(here::here("processed-data","Metadata.csv")) # -- this is actually all samples


# 1.1 Quick disease fixes ----

## fix metadata so diseases are clearer
metadata$`Disease Simplified`[metadata$`Disease Simplified` == "None"] <- "Control" #None is disease state (so healthy control) GSE232753
metadata$`Disease Simplified`[metadata$`Disease Simplified` == "Healthy"] <- "Control"
metadata$`Disease Simplified`[metadata$`Disease Simplified` == "Healthy control"] <- "Control"
metadata$`Disease Simplified`[metadata$`Disease Simplified` == "No Sepsis"] <- "Potential for sepsis but was not developed" #Potential Sepsis, stil were sick GSE63311
metadata$`Disease Simplified`[metadata$`Disease Simplified` == "No-Sepsis"] <- "Critical Illness" # Critical Illness from GSE189400
metadata$`Disease Simplified`[metadata$`Disease Simplified` == "Noninfectious"] <- "ICU noninfectious condition" #GSE134347
metadata$`Disease Simplified`[metadata$`Disease Simplified` == "Sepsis with ARDs"] <- "Sepsis - ARDs"
metadata$`Disease Simplified`[metadata$`Disease Simplified` == "Septic Shock with ARDs"] <- "Sepsis - ARDs"
metadata$`Disease Simplified`[metadata$`Disease Simplified` == "Uncomplicated Infection"] <- "Uncomplicated baterical infection" #GSE154918
metadata$`Disease Simplified`[metadata$`Disease Simplified` == "untreated"] <- "Intubated subjects undergoing mechanical ventilation" #GSE32707
metadata$`Disease Simplified`[metadata$`Disease Simplified` == "Septic ARDs"] <- "Sepsis - ARDs"
metadata$`Disease Simplified`[metadata$`Disease Simplified` == "Septic Shock"] <- "Sepsis - Shock"
metadata$`Disease Simplified`[metadata$`Disease Simplified` == "sepsis"] <- "Sepsis" 
metadata$`Disease Simplified`[metadata$`Disease Simplified` == "normal"] <- "Control" 
metadata$`Disease Simplified`[metadata$`Disease Simplified` == "SAKI"] <- "Sepsis - AKI" 
dim(metadata)
metadata <- metadata[!is.na(metadata$`Disease Simplified`),] #we lose 350 samples here with NA values

## Getting rid of these to be safe - not sure why they are excluded, not excluding any other samples for now
metadata <-  metadata[!(metadata$`Disease Simplified` == "Excluded"),] #lose 5 samples here
dim(metadata)

# 1.2 Create levels of all diseases ----
metadata$`Disease Simplified` <- factor(metadata$`Disease Simplified`,levels = c("Control","Intubated subjects undergoing mechanical ventilation",
                                                                                 "Uncomplicated baterical infection","Anaphylaxis",
                                                                                 "ICU noninfectious condition","Head Trauma",
                                                                                 "Critical Illness","COVID-19","Cardiogenic Shock",
                                                                                 "Potential for sepsis but was not developed","SIRS",
                                                                                 "Sepsis - AKI","Sepsis","Sepsis - ARDs","Sepsis - Shock"))


## * -- save RDS -- *
saveRDS(metadata, here::here("data","metadata.rds"))
