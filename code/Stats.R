

# 1.0 Read in Data ------------------------------------------------------------
source(here::here("final_code","ReadData.R"))


meta <- readDiseaseMeta()
head(meta)
dim(meta)
table(meta$Survival)
table(meta$Shock)
meta_cut <- meta |> dplyr::select(`Sample ID`, Dataset ,Survival, Race,Gender,Age,
                                  `Disease Simplified`, PreviouslyDefinedSubtypes,ClusterID)
head(meta_cut)

# meta_cut$`Disease Simplified`[meta_cut$`Disease Simplified` == "Sepsis - ARDs"] <- "Sepsis"
# meta_cut$`Disease Simplified`[meta_cut$`Disease Simplified` == "Sepsis - AKI"] <- "Sepsis"
# relevel disease
meta_cut$`Disease Simplified` <- factor(meta_cut$`Disease Simplified`, levels = c("Sepsis","Sepsis - Shock"))

## Without proportions: here::here("results","KM","Stats","ChiSquaredTests.tsv")

# Lets try to get the null expected distribution of each cluster given  --------
# These hopwfully will help tell us expected proportion to include missingness for all samples
meta_cut_alt <- meta_cut |> dplyr::rowwise() |>
  dplyr::mutate(SurviveNull = sample(c("Survived","Died"), 1))

meta_cut_alt <- meta_cut_alt |>dplyr::rowwise() |>
  dplyr::mutate(GenderNull = sample(c("Female","Male"), 1))

meta_cut_alt <- meta_cut_alt|> dplyr::rowwise() |>
  dplyr::mutate(RaceNull = sample(c("African American", "Asian","Caucasian", "Other"), 1))

# -- Endotypes --
meta_cut_alt <- meta_cut_alt |> dplyr::rowwise() |>
  dplyr::mutate(HyperNull = sample(c("Hyper","Hypo"), 1))

meta_cut_alt <- meta_cut_alt |> dplyr::rowwise() |>
  dplyr::mutate(LipidNull = sample(c("Lipid_Hypo","Lipid_Normo"), 1))

meta_cut_alt <- meta_cut_alt |> dplyr::rowwise() |>
  dplyr::mutate(SweeneyNull = sample(c("Inflammopathic", "Adaptive","Coagulopathic","Unclustered"), 1))

meta_cut_alt <- meta_cut_alt |> dplyr::rowwise() |>
  dplyr::mutate(MARSNull = sample(c("Mars1", "Mars2","Mars3","Mars4"), 1))

write.table(meta_cut_alt, here::here("final_files","MetadataFinalStatsNulls.tsv"), sep = '\t', col.names = TRUE,row.names = FALSE)


# 2.0 Chi-Squared Tests ----


# SURVIVAL ----------------------------------------------------------------

# First the null
sink(here::here("results","KM","Stats","Survival_ChiSquaredKM22.log"), append = TRUE)
cat(" ---- NULL DISTRIBUTION ----\n")
CR <- chisq.test(meta_cut_alt$SurviveNull, meta_cut_alt$ClusterID)
cat("-- expected distribution -- \n")
CR$expected
cat("-- observed distribution -- \n")
table(meta_cut_alt$SurviveNull, meta_cut_alt$ClusterID)
cat("\n-- X stat -- \n")
CR$statistic
cat("\n-- P value Before Bonferroni Correction:-- \n")
CR$p.value
cat("\n-- P value After Bonferroni Correction:-- \n")
p.adjust(chisq.test(meta_cut_alt$SurviveNull, meta_cut_alt$ClusterID)$p.value, method = "bonferroni", n = 10)
cat("# -- Proportion: -- #\n")
propsurv <- as.data.frame(table(meta_cut_alt$SurviveNull, meta_cut_alt$ClusterID))
propsurv |>
  dplyr::group_by(Var1) |>
  dplyr::mutate(Prop = Freq / sum(Freq)) |>
  dplyr::ungroup()

cat("\n\n ---- TRUE SAMPLES ----\n")
survive <- meta_cut_alt[meta_cut_alt$Survival %in% c("Died","Survived"),]
cat("## ---- Survival ---- ##\n")
cat("-- dim -- \n")
dim(survive)
CR <- chisq.test(survive$Survival, survive$ClusterID)
CR
cat("-- expected distribution -- \n")
CR$expected
cat("-- observed distribution -- \n")
table(survive$Survival, survive$ClusterID)
cat("\n-- X stat -- \n")
CR$statistic
cat("\n-- P value Before Bonferroni Correction:-- \n")
CR$p.value
cat("\n-- P value After Bonferroni Correction:-- \n")
p.adjust(chisq.test(survive$Survival, survive$ClusterID)$p.value, method = "bonferroni", n = 10)
cat("# -- Proportion: -- #\n")
propsurv <- as.data.frame(table(survive$Survival, survive$ClusterID))
propsurv |>
  dplyr::group_by(Var1) |>
  dplyr::mutate(Prop = Freq / sum(Freq)) |>
  dplyr::ungroup()

sink()


# 2.2 Gender ----

sink(here::here("results","KM","Stats","Gender_ChiSquaredKM22.log"), append = TRUE)
cat(" ---- NULL DISTRIBUTION ----\n")
CR <- chisq.test(meta_cut_alt$GenderNull, meta_cut_alt$ClusterID)
CR
cat("-- expected distribution -- \n")
CR$expected
cat("-- observed distribution -- \n")
table(meta_cut_alt$GenderNull, meta_cut_alt$ClusterID)
cat("\n-- X stat -- \n")
CR$statistic
cat("\n-- P value Before Bonferroni Correction:-- \n")
CR$p.value
cat("\n-- P value After Bonferroni Correction:-- \n")
p.adjust(chisq.test(meta_cut_alt$GenderNull, meta_cut_alt$ClusterID)$p.value, method = "bonferroni", n = 10)
cat("# -- Proportion: -- #\n")
propsurv <- as.data.frame(table(meta_cut_alt$GenderNull, meta_cut_alt$ClusterID))
propsurv |>
  dplyr::group_by(Var1) |>
  dplyr::mutate(Prop = Freq / sum(Freq)) |>
  dplyr::ungroup()

cat("\n\n ---- TRUE SAMPLES ----\n")
all_gender <- meta_cut_alt[meta_cut_alt$Gender %in% c("Female","Male"),]
cat("## ---- Gender ---- ##\n")
cat("-- dim -- \n")
dim(all_gender)
CR<-chisq.test(all_gender$Gender, all_gender$ClusterID)
CR
cat("-- expected distribution -- \n")
CR$expected
cat("-- observed distribution -- \n")
table(all_gender$Gender, all_gender$ClusterID)
cat("\n-- X stat -- \n")
CR$statistic
cat("\n-- P value Before Bonferroni Correction:-- \n")
CR$p.value
cat("\n-- P value After Bonferroni Correction:-- \n")
p.adjust(chisq.test(all_gender$Gender, all_gender$ClusterID)$p.value, method = "bonferroni", n = 10)
cat("# -- Proportion -- #")
propgender <- as.data.frame(table(all_gender$Gender, all_gender$ClusterID))
propgender |>
  dplyr::group_by(Var1) |>
  dplyr::mutate(Prop = Freq / sum(Freq)) |>
  dplyr::ungroup()
sink()




# AGE ---------------------------------------------------------------------
sink(here::here("results","KM","Stats","Age_ANOVAKM22.log"), append = TRUE)

age_cut <- meta[!(is.na(meta$Age)),] |> dplyr::select(Age, ClusterID, Dataset)
age_cut <- age_cut[!(age_cut$Age %in% c("20-39","> 80", "40-59", "60-79")),]
age_cut$Age <- as.numeric(age_cut$Age)
#get_age <- age_cut[age_cut$ClusterID]

cat("---- One-Way Anova ----")
oneway.test(Age ~ ClusterID,
            data = age_cut,
            var.equal = TRUE # assuming equal variances
)
age_test <- oneway.test(Age ~ ClusterID,
                        data = age_cut,
                        var.equal = TRUE # assuming equal variances
)

cat("---- Adjusted P-Value ----")
p.adjust(age_test$p.value, method = "bonferroni", n = 10)

sink()



# Disease ---------------------------------------------------------------------
sink(here::here("results","KM","Stats","Disease_ChiSquaredKM22.log"), append = TRUE)
# Disease does not need a null distribution

cat("\n\n ---- TRUE SAMPLES ----\n")
cat("## ---- Disease ---- ##\n")
cat("-- dim --")
dim(meta_cut)
CR <- chisq.test(meta_cut_alt$`Disease Simplified`, meta_cut_alt$ClusterID)
CR
cat("-- expected distribution -- \n")
CR$expected
cat("-- expected proportion -- \n")
exp_table <- t(CR$expected) |> as.data.frame()
exp_long  <- tidyr::pivot_longer(exp_table, cols = c("Sepsis","Sepsis - Shock"), names_to = "Disease",values_to = "ExpectedCount")
exp_long |> 
  dplyr::group_by(Disease) |>
  dplyr::mutate(Prop = ExpectedCount / sum(ExpectedCount)) |>
  dplyr::ungroup()
cat("-- observed distribution -- \n")
table(meta_cut_alt$`Disease Simplified`, meta_cut_alt$ClusterID)
cat("\n-- X stat -- \n")
CR$statistic
cat("\n-- P value Before Bonferroni Correction:-- \n")
CR$p.value
cat("\n-- P value After Bonferroni Correction:-- \n")
p.adjust(chisq.test(meta_cut_alt$`Disease Simplified`, meta_cut_alt$ClusterID)$p.value, method = "bonferroni", n = 10)
cat("# -- Proportion -- #")
propdisease <- as.data.frame(table(meta_cut_alt$`Disease Simplified`, meta_cut_alt$ClusterID))
propdisease |>
  dplyr::group_by(Var1) |>
  dplyr::mutate(Prop = Freq / sum(Freq)) |>
  dplyr::ungroup()
sink()


# RACE --------------------------------------------------------------------
sink(here::here("results","KM","Stats","Race_ChiSquaredKM22.log"), append = TRUE)
cat(" ---- NULL DISTRIBUTION ----\n")
CR <- chisq.test(meta_cut_alt$RaceNull, meta_cut_alt$ClusterID)
cat("-- expected distribution -- \n")
CR$expected
CR
cat("-- observed distribution -- \n")
table(meta_cut_alt$RaceNull, meta_cut_alt$ClusterID)
cat("\n-- X stat -- \n")
CR$statistic
cat("\n-- P value Before Bonferroni Correction:-- \n")
CR$p.value
cat("\n-- P value After Bonferroni Correction:-- \n")
p.adjust(chisq.test(meta_cut_alt$RaceNull, meta_cut_alt$ClusterID)$p.value, method = "bonferroni", n = 10)
cat("# -- Proportion: -- #\n")
propsurv <- as.data.frame(table(meta_cut_alt$RaceNull, meta_cut_alt$ClusterID))
propsurv |>
  dplyr::group_by(Var1) |>
  dplyr::mutate(Prop = Freq / sum(Freq)) |>
  dplyr::ungroup()

cat("\n\n ---- TRUE SAMPLES ----\n")
race <- meta_cut_alt[meta_cut_alt$Race %in% c("African American", "Asian","Caucasian", "Other"),]
cat("## ---- Race ---- ##\n")
cat("-- dim --")
dim(race) #313 had race info -

CR<-chisq.test(race$Race, race$ClusterID)
CR
cat("-- expected distribution -- \n")
CR$expected
cat("-- observed distribution -- \n")
table(race$Race, race$ClusterID) #  # TODO we need to replace this with Fishers
cat("\n-- X stat -- \n")
CR$statistic
cat("\n-- P value Before Bonferroni Correction:-- \n")
CR$p.value
cat("\n-- P value After Bonferroni Correction:-- \n")
p.adjust(chisq.test(race$Race, race$ClusterID)$p.value, method = "bonferroni", n = 10)
cat("# -- Proportion -- #")
proprace <- as.data.frame(table(race$Race, race$ClusterID))
proprace |>
  dplyr::group_by(Var1) |>
  dplyr::mutate(Prop = Freq / sum(Freq)) |>
  dplyr::ungroup()
sink()

## Fishers exact 
sink(here::here("results","KM","Stats","Race_ChiSquaredKM22.log"), append = TRUE)
# 
cat("## ---- Race FISHERS---- ##\n")
cat("-- dim --")
dim(race) #313 had race info -

#fisher.test(table(race$Race, race$ClusterID),workspace = 2e7)
CR<-fisher.test(table(race$Race, race$ClusterID),simulate.p.value=TRUE)
CR
cat("-- expected distribution -- \n")
CR$expectedj
cat("-- observed distribution -- \n")
table(race$Race, race$ClusterID) #  # TODO we need to replace this with Fishers
cat("\n-- X stat -- \n")
CR$statistic
cat("\n-- P value Before Bonferroni Correction:-- \n")
CR$p.value
cat("\n-- P value After Bonferroni Correction:-- \n")
p.adjust(fisher.test(table(race$Race, race$ClusterID),simulate.p.value=TRUE)$p.value, method = "bonferroni", n = 10)
cat("# -- Proportion -- #")
proprace <- as.data.frame(table(race$Race, race$ClusterID))
proprace |>
  dplyr::group_by(Var1) |>
  dplyr::mutate(Prop = Freq / sum(Freq)) |>
  dplyr::ungroup()
sink()



# 2.6 Endotypes: Hypo/Hyper ---- NEW BC N COUNTED

sink(here::here("results","KM","Stats","HyperHypo_ChiSquaredKM22.log"), append = TRUE)
cat(" ---- NULL DISTRIBUTION ----\n")
CR <- chisq.test(meta_cut_alt$HyperNull, meta_cut_alt$ClusterID)
CR
cat("-- expected distribution -- \n")
CR$expected
cat("-- observed distribution -- \n")
table(meta_cut_alt$HyperNull, meta_cut_alt$ClusterID)
cat("\n-- X stat -- \n")
CR$statistic
cat("\n-- P value Before Bonferroni Correction:-- \n")
CR$p.value
cat("\n-- P value After Bonferroni Correction:-- \n")
p.adjust(chisq.test(meta_cut_alt$HyperNull, meta_cut_alt$ClusterID)$p.value, method = "bonferroni", n = 10)
cat("# -- Proportion: -- #\n")
propsurv <- as.data.frame(table(meta_cut_alt$HyperNull, meta_cut_alt$ClusterID))
propsurv |>
  dplyr::group_by(Var1) |>
  dplyr::mutate(Prop = Freq / sum(Freq)) |>
  dplyr::ungroup()

cat("\n\n ---- TRUE SAMPLES ----\n")
hypo <- meta_cut_alt[meta_cut_alt$`PreviouslyDefinedSubtypes` %in% c("Hyper","Hypo"),]

cat("## ---- Endotypes: Hypo/Hyper ---- ##\n")
cat("-- dim -- \n")
dim(hypo)
CR<-chisq.test(hypo$`PreviouslyDefinedSubtypes`, hypo$ClusterID) # TODO we need to replace this with Fishers
CR
cat("-- expected distribution -- \n")
CR$expected
cat("-- observed distribution -- \n")
table(hypo$`PreviouslyDefinedSubtypes`, hypo$ClusterID)
cat("\n-- X stat -- \n")
CR$statistic
cat("\n-- P value Before Bonferroni Correction:-- \n")
CR$p.value
cat("\n-- P value After Bonferroni Correction:-- \n")
p.adjust(chisq.test(hypo$`PreviouslyDefinedSubtypes`, hypo$ClusterID)$p.value, method = "bonferroni", n = 10)
cat("# -- Proportion -- #")
propendo <- as.data.frame(table(hypo$`PreviouslyDefinedSubtypes`, hypo$ClusterID))
propendo |>
  dplyr::group_by(Var1) |>
  dplyr::mutate(Prop = Freq / sum(Freq)) |>
  dplyr::ungroup()
sink()

# 2.7 Endotypes: Lipid Hypo/Normo ----

sink(here::here("results","KM","Stats","LipidHypoerNormo_ChiSquaredKM22.log"), append = TRUE)
cat(" ---- NULL DISTRIBUTION ----\n")
CR <- chisq.test(meta_cut_alt$LipidNull, meta_cut_alt$ClusterID)
cat("-- expected distribution -- \n")
CR$expected
cat("-- observed distribution -- \n")
table(meta_cut_alt$LipidNull, meta_cut_alt$ClusterID)
cat("\n-- X stat -- \n")
CR$statistic
cat("\n-- P value Before Bonferroni Correction:-- \n")
CR$p.value
cat("\n-- P value After Bonferroni Correction:-- \n")
p.adjust(chisq.test(meta_cut_alt$LipidNull, meta_cut_alt$ClusterID)$p.value, method = "bonferroni", n = 10)
cat("# -- Proportion: -- #\n")
propsurv <- as.data.frame(table(meta_cut_alt$LipidNull, meta_cut_alt$ClusterID))
propsurv |>
  dplyr::group_by(Var1) |>
  dplyr::mutate(Prop = Freq / sum(Freq)) |>
  dplyr::ungroup()

cat("\n\n ---- TRUE SAMPLES ----\n")
Lipid <- meta_cut_alt[meta_cut_alt$`PreviouslyDefinedSubtypes` %in% c("Lipid_Normo","Lipid_Hypo"),]

cat("## ---- Endotypes: Lipid Hypo/Normo ---- ##\n")
cat("-- dim --")
dim(Lipid)
CR<-chisq.test(Lipid$`PreviouslyDefinedSubtypes`, Lipid$ClusterID)
CR
cat("-- expected distribution -- \n")
CR$expected
cat("-- observed distribution -- \n")
table(Lipid$`PreviouslyDefinedSubtypes`, Lipid$ClusterID)
cat("\n-- X stat -- \n")
CR$statistic
cat("\n-- P value Before Bonferroni Correction:-- \n")
CR$p.value
cat("\n-- P value After Bonferroni Correction:-- \n")
p.adjust(chisq.test(Lipid$`PreviouslyDefinedSubtypes`, Lipid$ClusterID)$p.value, method = "bonferroni", n = 10)
cat("# -- Proportion -- #")
propLipid <- as.data.frame(table(Lipid$`PreviouslyDefinedSubtypes`, Lipid$ClusterID))
propLipid |>
  dplyr::group_by(Var1) |>
  dplyr::mutate(Prop = Freq / sum(Freq)) |>
  dplyr::ungroup()
sink()

## Fishers exact 
sink(here::here("results","KM","Stats","LipidHypoerNormo_ChiSquaredKM22.log"), append = TRUE)
# 
cat("## ---- Lipid Endotype FISHERS---- ##\n")
cat("-- dim --")
dim(Lipid)
CR<-fisher.test(Lipid$`PreviouslyDefinedSubtypes`, Lipid$ClusterID, simulate.p.value = TRUE)
CR
cat("-- expected distribution -- \n")
CR$expected
cat("-- observed distribution -- \n")
table(Lipid$`PreviouslyDefinedSubtypes`, Lipid$ClusterID)
cat("\n-- X stat -- \n")
CR$statistic
cat("\n-- P value Before Bonferroni Correction:-- \n")
CR$p.value
cat("\n-- P value After Bonferroni Correction:-- \n")
p.adjust(fisher.test(Lipid$`PreviouslyDefinedSubtypes`, Lipid$ClusterID, simulate.p.value = FALSE)$p.value, method = "bonferroni", n = 10)
cat("# -- Proportion -- #")
propLipid <- as.data.frame(table(Lipid$`PreviouslyDefinedSubtypes`, Lipid$ClusterID))
propLipid |>
  dplyr::group_by(Var1) |>
  dplyr::mutate(Prop = Freq / sum(Freq)) |>
  dplyr::ungroup()
sink()



# 2.8 Endotypes: Mars ----


sink(here::here("results","KM","Stats","MARs_ChiSquaredKM22.log"), append = TRUE)
cat(" ---- NULL DISTRIBUTION ----\n")
CR <- chisq.test(meta_cut_alt$MARSNull, meta_cut_alt$ClusterID)
cat("-- expected distribution -- \n")
CR
CR$expected
cat("-- observed distribution -- \n")
table(meta_cut_alt$MARSNull, meta_cut_alt$ClusterID)
cat("\n-- X stat -- \n")
CR$statistic
cat("\n-- P value Before Bonferroni Correction:-- \n")
CR$p.value
cat("\n-- P value After Bonferroni Correction:-- \n")
p.adjust(chisq.test(meta_cut_alt$MARSNull, meta_cut_alt$ClusterID)$p.value, method = "bonferroni", n = 10)
cat("# -- Proportion: -- #\n")
propsurv <- as.data.frame(table(meta_cut_alt$MARSNull, meta_cut_alt$ClusterID))
propsurv |>
  dplyr::group_by(Var1) |>
  dplyr::mutate(Prop = Freq / sum(Freq)) |>
  dplyr::ungroup()

cat("\n\n ---- TRUE SAMPLES ----\n")
mars <- meta_cut[grep("Mars",meta_cut$`PreviouslyDefinedSubtypes`),]
cat("## ---- Endotypes: Mars ---- ##\n")
cat("-- dim --")
dim(mars)
CR<-chisq.test(mars$`PreviouslyDefinedSubtypes`, mars$ClusterID)
CR
cat("-- expected distribution -- \n")
CR$expected
cat("-- observed distribution -- \n")
table(mars$`PreviouslyDefinedSubtypes`, mars$ClusterID)
cat("\n-- X stat -- \n")
CR$statistic
cat("\n-- P value Before Bonferroni Correction:-- \n")
CR$p.value
cat("\n-- P value After Bonferroni Correction:-- \n")
p.adjust(chisq.test(mars$`PreviouslyDefinedSubtypes`, mars$ClusterID)$p.value, method = "bonferroni", n = 10)
cat("# -- Proportion -- #")
propmars <- as.data.frame(table(mars$`PreviouslyDefinedSubtypes`, mars$ClusterID))
propmars |>
  dplyr::group_by(Var1) |>
  dplyr::mutate(Prop = Freq / sum(Freq)) |>
  dplyr::ungroup()
sink()


# 2.9 Endotype: Sweeney
sink(here::here("results","KM","Stats","Sweeney_ChiSquaredKM22.log"), append = TRUE)
cat(" ---- NULL DISTRIBUTION ----\n")
CR<- chisq.test(meta_cut_alt$SweeneyNull, meta_cut_alt$ClusterID)
CR
cat("-- expected distribution -- \n")
CR$expected
cat("-- observed distribution -- \n")
table(meta_cut_alt$SweeneyNull, meta_cut_alt$ClusterID)
cat("\n-- X stat -- \n")
CR$statistic
cat("\n-- P value Before Bonferroni Correction:-- \n")
CR$p.value
cat("\n-- P value After Bonferroni Correction:-- \n")
p.adjust(chisq.test(meta_cut_alt$SweeneyNull, meta_cut_alt$ClusterID)$p.value, method = "bonferroni", n = 10)
cat("# -- Proportion: -- #\n")
propsurv <- as.data.frame(table(meta_cut_alt$SweeneyNull, meta_cut_alt$ClusterID))
propsurv |>
  dplyr::group_by(Var1) |>
  dplyr::mutate(Prop = Freq / sum(Freq)) |>
  dplyr::ungroup()

cat("\n\n ---- TRUE SAMPLES ----\n")
cat("## ---- Endotype: Sweeney ---- ##\n")
cat("-- dim --")
sweeney <-  meta_cut_alt[meta_cut_alt$`PreviouslyDefinedSubtypes` %in% c("Adaptive","Coagulopathic","Inflammopathic","Unclustered"),]
CR<-chisq.test(sweeney$PreviouslyDefinedSubtypes, sweeney$ClusterID)
CR
cat("-- expected distribution -- \n")
CR$expected
cat("-- observed distribution -- \n")
table(sweeney$PreviouslyDefinedSubtypes, sweeney$ClusterID) # TODO we need to replace this with Fishers
cat("\n-- X stat -- \n")
CR$statistic
cat("\n-- P value Before Bonferroni Correction:-- \n")
CR$p.value
cat("\n-- P value After Bonferroni Correction:-- \n")
p.adjust(chisq.test(sweeney$PreviouslyDefinedSubtypes, sweeney$ClusterID)$p.value, method = "bonferroni", n = 10)
cat("# -- Proportion -- #")
proprace <- as.data.frame(table(sweeney$PreviouslyDefinedSubtypes, sweeney$ClusterID))
proprace |>
  dplyr::group_by(Var1) |>
  dplyr::mutate(Prop = Freq / sum(Freq)) |>
  dplyr::ungroup()
sink()

sink(here::here("results","KM","Stats","Sweeney_ChiSquaredKM22.log"), append = TRUE)
cat("## ---- Endotype: Sweeney FISHERS ---- ##\n")
cat("-- dim --")
dim(sweeney) #337

CR<-fisher.test(sweeney$PreviouslyDefinedSubtypes, sweeney$ClusterID, simulate.p.value = TRUE)
CR
cat("-- expected distribution -- \n")
CR$expected
cat("-- observed distribution -- \n")
table(sweeney$PreviouslyDefinedSubtypes, sweeney$ClusterID) # TODO we need to replace this with Fishers
cat("\n-- X stat -- \n")
CR$statistic
cat("\n-- P value Before Bonferroni Correction:-- \n")
CR$p.value
cat("\n-- P value After Bonferroni Correction:-- \n")
p.adjust(fisher.test(sweeney$PreviouslyDefinedSubtypes, sweeney$ClusterID, simulate.p.value = TRUE)$p.value, method = "bonferroni", n = 10)
cat("# -- Proportion -- #")
proprace <- as.data.frame(table(sweeney$PreviouslyDefinedSubtypes, sweeney$ClusterID))
proprace |>
  dplyr::group_by(Var1) |>
  dplyr::mutate(Prop = Freq / sum(Freq)) |>
  dplyr::ungroup()
sink()



# All Endotytpes for clinical barchaart -----------------------------------

sink(here::here("results","KM","Stats","AllEndotypes_ChiSquaredKM22.log"), append = TRUE)

cat("\n\n ---- TRUE SAMPLES ----\n")
all <- meta_cut[!(is.na(meta_cut$PreviouslyDefinedSubtypes)),]
all$Study <- NA
all$Study[all$PreviouslyDefinedSubtypes %in% c("Coagulopathic","Unclustered","Inflammopathic","Adaptive")] <- "Sweeney"
all$Study[all$PreviouslyDefinedSubtypes %in% c("Mars1","Mars2","Mars3","Mars4")] <- "Mars"
all$Study[all$PreviouslyDefinedSubtypes %in% c("Lipid_Hypo","Lipid_Normo")] <- "Guirgis"
all$Study[all$PreviouslyDefinedSubtypes %in% c("Hypo","Hyper")] <- "Calfee"

cat("## ---- Endotypes: all chi squared ---- ##\n")
cat("-- dim --")
dim(all)
CR<-chisq.test(all$Study, all$ClusterID)
CR
cat("-- expected distribution -- \n")
CR$expected
cat("-- observed distribution -- \n")
table(all$Study, all$ClusterID)
cat("\n-- X stat -- \n")
CR$statistic
cat("\n-- P value Before Bonferroni Correction:-- \n")
CR$p.value
cat("\n-- P value After Bonferroni Correction:-- \n")
p.adjust(chisq.test(all$Study, all$ClusterID)$p.value, method = "bonferroni", n = 10)
cat("# -- Proportion -- #")
propall <- as.data.frame(table(all$Study, all$ClusterID))
propall |>
  dplyr::group_by(Var1) |>
  dplyr::mutate(Prop = Freq / sum(Freq)) |>
  dplyr::ungroup()
# 
# cat("---- FISHERS ----")
# p.adjust(fisher.test(all$PreviouslyDefinedSubtypes,all$ClusterID, simulate.p.value = TRUE)$p.value, method = "bonferroni", n = 10)
# cat("# -- Proportion -- #")
# proprace <- as.data.frame(table(all$PreviouslyDefinedSubtypes, all$ClusterID))
# proprace |>
#   dplyr::group_by(Var1) |>
#   dplyr::mutate(Prop = Freq / sum(Freq)) |>
#   dplyr::ungroup()
sink()
