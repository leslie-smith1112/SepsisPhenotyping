## Perform batch correction on all our studies.


# 1.0 Read in Expression and Metadata ---- 

## Read in batchIDs, metadata and expression
batchIDs <- readRDS(here::here("data","batchIDs_allSamples.rds"))
metadata <- readRDS(here::here("data","all_samples_metadata.rds"))
all_expr <- readRDS(here::here("data","data","raw_all_samples_expression.rds"))

## check dimensions all_expr
dim(all_expr)
all_expr[1:5,1:5]


# 2.0 Ensure Same Samples in all Data ----

## make sure batchIDs, metadata, and expression matrix all have same samples 
## - keep metadata only from samples we have expression data for - #
dim(metadata)#5108 -> excluding NA 4753
meta.data <- metadata[metadata$`Sample ID` %in% colnames(all_expr),]
dim(meta.data)#4753 samples -> 3754 
dim(all_expr)
all_expr <- all_expr[,colnames(all_expr) %in% meta.data$`Sample ID`]
dim(all_expr) #3754
all_expr[1:5,1:5]
batchIDs <- batchIDs[names(batchIDs) %in% colnames(all_expr)] # I think we lose that sample Dongyuan mentioned here
length(batchIDs)

## Double checking
all(meta.data$`Sample ID` %in% colnames(all_expr)) 
all(colnames(all_expr) %in% meta.data$`Sample ID` ) 
## the order of these must be equal:
all.equal(colnames(all_expr), names(batchIDs)) 


# 3.0 Plot Initial PCA of All Samples ----

## Trasnpose matrix for PCA function (want PCA of samples not genes)
all_transposed <- t(all_expr) 
all_transposed[1:5,1:5]
dim(all_transposed)

## - get rid of 0 values
all_pcaMat <- all_transposed[ , which(apply(all_transposed, 2, var) != 0)]


## 3.1 PCA Call: ----

all_pca <- prcomp(all_pcaMat, center = TRUE, scale. = TRUE)
dim(all_pca$x) #

## * -- save RDS -- *
saveRDS(all_pca$x, here::here("data","Init_PCA_values_allSamples.rds"))


## 3.2  Plot results ----
library(ggplot2)

#create matrix 
init_dtp <- data.frame("SampleID" = meta.data$`Sample ID`,"Series"=meta.data$Dataset, 'Disease' = meta.data$`Disease Simplified`, "Technology" = meta.data$`RNA-seq/Microarray`)
dat_pca_init <- tibble::rownames_to_column(as.data.frame(all_pca$x), "SampleID")
dtp <- merge(init_dtp, dat_pca_init, by = "SampleID")

# -- Plotting -- 
p <- ggplot(data=dtp, aes(x=PC1, y=PC2,  color=Series,shape = Technology)) + geom_point() + 
  ggtitle(paste0("All ",nrow(meta.data)," Samples - PCA No Batch Correction by Series")) + 
  labs(y= "PC2", x = "PC1")
p
ggsave(here::here("figures","Initial_PCA_Dataset.png"),plot = last_plot(),width = 7, height = 7)

# -- Plotting --
p <- ggplot(data=dtp, aes(x=PC1, y=PC2,  color=Disease,shape = Technology)) + geom_point(aes(color = Disease, alpha= 0.8)) + 
  ggtitle(paste0("All ",nrow(meta.data)," Samples, ", nrow(all), "Genes - PCA No Batch Correction by Disease")) + 
  labs(y= "PC2", x = "PC1") + 
  scale_color_manual(values = c("#f20909", "#8f8d98","#b6b5b6","#a48793","#6a608f","#576d78",
                               "#572d57", "#89063b","#4b795f","#572111","#7f8061","#07d81f", "#36bfc7","#3e85d0","#573ed0")) +
  theme_bw()
p
ggsave(here::here("figures","Initial_PCA_Disease.png"),plot = last_plot(),width = 8, height = 7)

 
# 4.0 Batch Correction  with All Samples ----

## Double checking dimensions and diseases
all_expr[1:5,1:5]
dim(meta.data)
dim(all_expr)
table(meta.data$`Disease Simplified`)
all.equal(colnames(all_expr), names(batchIDs))


## 4.1 Create Model Matrix ----

design <- model.matrix(~`Disease Simplified`, meta.data)
head(design)


## 4.2 Create Model Matrix ----


## get rid of negative values for ComBat 
raw_merged <- as.matrix(all_expr - min(all_expr))


## 4.3 Run ComBat ----
library(sva)
dat_batch_adjusted_norm_new <- ComBat(raw_merged, batchIDs, mod = design)
dat_batch_adjusted_norm_new[1:5,1:5]

## * -- save RDS -- *
saveRDS(dat_batch_adjusted_norm_new, here::here("data","all_samples_expression.rds"))



## trasnpose matrix for PCA function 
dat_adjusted_norm_transposed <- t(dat_batch_adjusted_norm_new)


## 4.4 Create PCA of Batch Corrected Samples

# get rid of 0 values 
dat_adjusted_norm_pcaMat <- dat_adjusted_norm_transposed[ , which(apply(dat_adjusted_norm_transposed, 2, var) != 0)]
dat_pca <- prcomp(dat_adjusted_norm_pcaMat, center = TRUE, scale. = TRUE)
dat_pca$x[1:5,1:5]

## * -- save RDS -- * 
saveRDS(dat_pca$x, here::here("data","BatchCorrected_PCA_values_allSamples.rds"))


## 4.5 Plot Results ----

## 4.5.1 Plot All Samples ----

# create matrix with PCA
init_dtp <- data.frame("SampleID" = meta.data$`Sample ID`,"Series"=meta.data$Dataset, 'Disease' = meta.data$`Disease Simplified`, "Technology" = meta.data$`RNA-seq/Microarray`)
dat_pca_init <- tibble::rownames_to_column(as.data.frame(dat_pca$x), "SampleID")
dtp <- merge(init_dtp, dat_pca_init, by = "SampleID")

# -- Plotting -- 
q <- ggplot(data=dtp, aes(x=PC1, y=PC2, color=Series, shape = Technology)) + geom_point() + ggtitle(paste0("Sepsis PCA Batch Corrected All Samples (",ncol(dat_batch_adjusted_norm_new),") Nonnorm Series")) + labs(y= "PC2", x = "PC1")
q
ggsave(here::here("figures","BatchCorrected_PCA_Dataset.png"),plot = last_plot(),width = 7, height = 7)

# -- Plotting -- 
p <- ggplot(data=dtp, aes(x=PC1, y=PC2,  color=Disease,shape = Technology)) + geom_point(aes(color = Disease, alpha= 0.8)) + 
  ggtitle(paste0("All ",nrow(meta.data)," Samples - PCA Batch Correction by Disease")) + 
  labs(y= "PC2", x = "PC1") + 
  scale_color_manual(values = c("#f20909", "#8f8d98","#b6b5b6","#a48793","#6a608f","#576d78",
                                "#572d57", "#89063b","#4b795f","#572111","#7f8061","#07d81f", "#36bfc7","#3e85d0","#573ed0")) +
  theme_bw()
p
ggsave(here::here("figures","BatchCorrected_PCA_Disease.png"),plot = last_plot(),width = 10, height = 7)


## You have to run Metadata_Clean_PostBatchCorrection.R to get metadata_clean.rds before running this - didnt want to make an additional script ;/

## 4.5.2 Plot Disease Samples Only ----
meta_clean <- readRDS(here::here("data","disease_metadata.rds")) #this metadata has been fixed for uniformity
dim(meta_clean)
head(meta_clean)
expr_dat <- dat_batch_adjusted_norm_new[,colnames(dat_batch_adjusted_norm_new) %in% meta_clean$`Sample ID`]
dim(expr_dat)
expr_dat[1:5,1:5]

# * -- save RDS -- *
saveRDS(expr_dat, here::here("data","disease_expression_all_genes.rds"))

dtp_disease <- dtp[dtp$SampleID %in% meta_clean$`Sample ID`,]
dim(dtp_disease)
#levels(dtp$Disease)
#have to reset dtp_disease levels
dtp_disease$Disease <- factor(dtp_disease$Disease,levels = c("Sepsis", "Sepsis - Shock", "Sepsis - ARDs", "Sepsis - AKI" ))

# -- Plotting -- 
p <- ggplot(data=dtp_disease, aes(x=PC1, y=PC2,  color=Disease,shape = Technology)) + geom_point(aes(color = Disease, alpha= 0.8)) + 
  ggtitle(paste0("All ",nrow(dtp_disease)," Samples, - PCA Batch Correction by Disease")) + 
  labs(y= "PC2", x = "PC1") + 
  scale_color_manual(values = c("#e5e1cc", "#fa8619","#192dfa","#07d81f")) +
  theme_bw()
p
ggsave(here::here("figures","BatchCorrected_PCA_DiseaseOnly.png"),plot = last_plot(),width = 10, height = 7)











