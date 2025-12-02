# Look at PCAs with only immune genes

# 1.0 GO terms downloaded with the following constraints: ----

# -------- 
# Downloaded immune system genes from GO term for immune response GO:0006955. Homoe Sapiens, proteins. url(https://www.ebi.ac.uk/QuickGO/annotations?geneProductType=protein&taxonId=9606&taxonUsage=descendants&goId=GO:0006955&goUsage=descendants&goUsageRelationships=is_a,part_of,occurs_in,regulates)
# * 45,671 terms
# * In GO terms:
#   * for Select how you want to use the terms that you have chosen: include child terms
# * for Select which relationship an annotated term should have to the terms above: is_a, part_of, occurs_in, regulates
# -------

##  -- To compare the GO terms with the top 5000 most variable - we have to PCA top 5000, GO terms, and we'll also separate do lipids as well
### -- then we will add all together and look at the PCA


## 1.1 Read in expression data

expr_dat <- readRDS(here::here("data","disease_expression_all_genes.rds"))
expr_dat[1:5,1:5]
dim(expr_dat)

## 1.2 Read in metadata ----
meta_clean <- readRDS(here::here("data","disease_metadata.rds"))


# 2.0 Top 5000 most variable genes ----

## 2.1 Get the top 5000 most variable genes using mean absolute deviation ----

gene_variability <- apply(expr_dat, 1, mad)
ranked_variances <- sort(gene_variability, decreasing = TRUE)
head(ranked_variances)
top_genes <- names(ranked_variances)[1:5000]
head(top_genes)

# always checking dimensions before and after cutting something
dim(expr_dat)
# keep 5000 most variable genes in the expression matrix. 
expr_dat_var <- expr_dat[rownames(expr_dat) %in% top_genes,] 
expr_dat_var[1:5,1:5]
dim(expr_dat_var)

## * -- save RDS -- * 
saveRDS(expr_dat_var, here::here("data","disease_expression.rds"))


## 2.2 PCA  ----

expr_var_t <- t(expr_dat_var)
expr_var_norm <- expr_var_t[ , which(apply(expr_var_t, 2, var) != 0)]
expr_var_norm[1:5,1:5]
dat_pca_var <- prcomp(expr_var_norm, center = TRUE, scale. = TRUE)
dat_pca_var$x[1:5,1:5]


# Create data frame for plotting
init_dtp <- data.frame("SampleID" = meta_clean$`Sample ID`,"Series"=meta_clean$Dataset, 'Disease' = meta_clean$`Disease Simplified`, "Technology" = meta_clean$`RNA-seq/Microarray`)
dat_pca_init <- tibble::rownames_to_column(as.data.frame(dat_pca_var$x), "SampleID")
dtp_var <- merge(init_dtp, dat_pca_init, by = "SampleID")
dtp_var[1:5,1:5]

## * -- save RDS -- * 
saveRDS(dtp_var, here::here("data","PCAValues_top5000Genes.rds"))


# -- Plotting --
library(ggplot2)
q <- ggplot(data=dtp_var, aes(x=PC1, y=PC2, color=Disease, shape = Technology)) + 
  ggtitle(paste0("5000 Most Variable Genes with ",nrow(dtp_var)," Samples")) + 
  geom_point(aes(alpha = 0.7)) + 
  scale_color_manual(values = c("#A59DA2","#a50a6c","#192dfa","#07d81f")) +
  theme_bw() + 
  labs(y= "PC2", x = "PC1")
q
ggsave(here::here("figures","PCAValues_5000MostVar.png"),plot = last_plot(),width = 7, height = 7)




# 3.0 Immune genes ----

gb_database <- readr::read_tsv(here::here("GO_Terms", "QuickGO-annotations.tsv"))
gb_database[1:5,1:5]
length(unique(gb_database$SYMBOL))
length(intersect(rownames(expr_dat), gb_database$SYMBOL)) #1207 intersecting genes
expr_dat[1:5,1:5]


## 3.1 Cut Genes to Immune Only ----
expr_imm <- expr_dat[rownames(expr_dat) %in% gb_database$SYMBOL,]
dim(expr_imm)
expr_imm[1:5,1:5]

## * -- save RDS -- * 
saveRDS(expr_imm, here::here("data","expr_immuneGenes.rds"))
dim(expr_imm)

# 3.2 PCA  ----

expr_imm_t <- t(expr_imm)
expr_imm_norm <- expr_imm_t[ , which(apply(expr_imm_t, 2, var) != 0)]
expr_imm_norm[1:5,1:5]
dat_pca_imm <- prcomp(expr_imm_norm, center = TRUE, scale. = TRUE)
dat_pca_imm$x[1:5,1:5]


# Create data frame for plotting
init_dtp <- data.frame("SampleID" = meta_clean$`Sample ID`,"Series"=meta_clean$Dataset, 'Disease' = meta_clean$`Disease Simplified`, "Technology" = meta_clean$`RNA-seq/Microarray`)
dat_pca_init <- tibble::rownames_to_column(as.data.frame(dat_pca_imm$x), "SampleID")
dtp_imm <- merge(init_dtp, dat_pca_init, by = "SampleID")
dtp_imm[1:5,1:5]

## * -- save RDS -- * 
saveRDS(dtp_imm, here::here("data","PCAValues_immuneGenes.rds"))


# -- Plotting --
q <- ggplot(data=dtp_imm, aes(x=PC1, y=PC2, color=Disease, shape = Technology)) + 
  geom_point(aes(alpha = 0.7)) + 
  scale_color_manual(values = c("#A59DA2","#a50a6c","#192dfa","#07d81f")) +
  theme_bw() + 
  ggtitle(paste0("GO Immune Genes with ",nrow(dtp_imm)," Samples")) + 
  labs(y= "PC2", x = "PC1")
q
ggsave(here::here("figures","PCAValues_immuneGenes.png"),plot = last_plot(),width = 7, height = 7)

# -- Plotting --
q <- ggplot(data=dtp_imm, aes(x=PC3, y=PC4, color=Disease, shape = Technology)) + 
  geom_point(aes(alpha = 0.7)) + 
  scale_color_manual(values = c("#A59DA2","#a50a6c","#192dfa","#07d81f")) +
  theme_bw() + 
  ggtitle(paste0("GO Immune Genes with ",nrow(dtp_imm)," Samples")) + 
  labs(y= "PC4", x = "PC3")
q
ggsave(here::here("figures","PCA34Values_immuneGenes.png"),plot = last_plot(),width = 7, height = 7)

# -- Plotting --
q <- ggplot(data=dtp_imm, aes(x=PC5, y=PC6, color=Disease, shape = Technology)) + 
  geom_point(aes(alpha = 0.7)) + 
  scale_color_manual(values = c("#A59DA2","#a50a6c","#192dfa","#07d81f")) +
  theme_bw() + 
  ggtitle(paste0("GO Immune Genes with ",nrow(dtp_imm)," Samples")) + 
  labs(y= "PC6", x = "PC5")
q
ggsave(here::here("figures","PCA56Values_immuneGenes.png"),plot = last_plot(),width = 7, height = 7)



# 4.0 Lipid genes ----

#Faheems lipid genes
genes <- c("CYP51A1","CYP4F2","LDLR","DHCR24","DHCR7","MSMO1","HMGCR","PCSK9","FDFT1","SQLE","ALOX5", "LSS","LCAT","LBR","CYP4Z1",
           "KCNH7", "CYP4X1","CYP46A1","APOA1","EBP","CYP4V2", "CYP39A1","SCARB1","HSD17B7","CYP4F22","ALOX15","ABCG1","NSDHL", "CYP4F12","ALOX12",
           "ABCA1","CYP4F11", "LOX","LIPE","TM7SF2","CYP4F8","PTGS1","LPA","LBR","CYP4F3","APOB","LBP","CETP","CYP4B1","CYP4A11",
           "FDFT1","LCAT","CYP4A22","PLTP","PTGS2")

## 4.1 Cut Genes to Immune Only ----
expr_lipid <- expr_dat[rownames(expr_dat) %in% genes,]
dim(expr_lipid)
expr_lipid[1:5,1:5]

## * -- save RDS -- * 
saveRDS(expr_lipid, here::here("data","expr_lipidGenes.rds"))


# 4.2 PCA  ----

expr_lipid_t <- t(expr_lipid)
expr_lipid_norm <- expr_lipid_t[ , which(apply(expr_lipid_t, 2, var) != 0)]
expr_lipid_norm[1:5,1:5]
dat_pca_lipid <- prcomp(expr_lipid_norm, center = TRUE, scale. = TRUE)
dat_pca_lipid$x[1:5,1:5]


# Create data frame for plotting
init_dtp <- data.frame("SampleID" = meta_clean$`Sample ID`,"Series"=meta_clean$Dataset, 'Disease' = meta_clean$`Disease Simplified`, "Technology" = meta_clean$`RNA-seq/Microarray`)
dat_pca_init <- tibble::rownames_to_column(as.data.frame(dat_pca_lipid$x), "SampleID")
dtp_lipid <- merge(init_dtp, dat_pca_init, by = "SampleID")
dtp_lipid[1:5,1:5]
dim(dtp_lipid)

## * -- save RDS -- * 
saveRDS(dtp_lipid, here::here("data","PCAValues_lipidGenes.rds"))


# -- Plotting --
q <- ggplot(data=dtp_lipid, aes(x=PC1, y=PC2, color=Disease, shape = Technology)) + 
  geom_point(aes(alpha = 0.7)) + 
  scale_color_manual(values = c("#A59DA2","#a50a6c","#192dfa","#07d81f")) +
  theme_bw() +  
  ggtitle(paste0("GO Lipid Genes with ",nrow(dtp_lipid)," Samples")) + 
  labs(y= "PC2", x = "PC1")
q
ggsave(here::here("figures","PCAValues_lipidGenes.png"),plot = last_plot(),width = 7, height = 7)



# 5.0 Create matrix with Immune - Lipid - and the rest most variable genes ---

go_genes <- unique(gb_database$SYMBOL)
length(go_genes)
keep_genes <- unique(c(genes, go_genes))
length(keep_genes)

#make sure genes are actually in the expression matrix that we have
keep_genes <- keep_genes[keep_genes %in% rownames(expr_dat)]
length(keep_genes) #1242
# use previous gene variance rankings
head(ranked_variances)
# variation selection should not have any overlapping genes with keep_genes
valid_var <- ranked_variances[!(names(ranked_variances) %in% keep_genes)]
head(valid_var)

#get the 5000 genes we will keep
all_genes <- c(keep_genes, names(valid_var[1:(5000 - length(keep_genes))]))
length(all_genes) #should both be 5000
length(unique(all_genes)) 

expr_mixed <- expr_dat[rownames(expr_dat) %in% all_genes,]
dim(expr_mixed)

## * -- save RDS -- * 
saveRDS(expr_mixed, here::here("data","expr_mixedGenes.rds"))


# 5.2 PCA  ----

expr_mixed_t <- t(expr_mixed)
expr_mixed_norm <- expr_mixed_t[ , which(apply(expr_mixed_t, 2, var) != 0)]
expr_mixed_norm[1:5,1:5]
dat_pca_mixed <- prcomp(expr_mixed_norm, center = TRUE, scale. = TRUE)
dat_pca_mixed$x[1:5,1:5]


# Create data frame for plotting
init_dtp <- data.frame("SampleID" = meta_clean$`Sample ID`,"Series"=meta_clean$Dataset, 'Disease' = meta_clean$`Disease Simplified`, "Technology" = meta_clean$`RNA-seq/Microarray`)
dat_pca_init <- tibble::rownames_to_column(as.data.frame(dat_pca_mixed$x), "SampleID")
dtp_mixed <- merge(init_dtp, dat_pca_init, by = "SampleID")
dtp_mixed[1:5,1:5]
dim(dtp_mixed)

## * -- save RDS -- * 
saveRDS(dtp_mixed, here::here("data","PCAValues_mixedGenes.rds"))


# -- Plotting --
q <- ggplot(data=dtp_mixed, aes(x=PC1, y=PC2, color=Disease, shape = Technology)) + 
  geom_point(aes(alpha = 0.7)) + 
  scale_color_manual(values = c("#A59DA2","#a50a6c","#192dfa","#07d81f")) +
  theme_bw() +  
  ggtitle(paste0("GO + Lipid + Variance Genes with ",nrow(dtp_mixed)," Samples")) + 
  labs(y= "PC2", x = "PC1")
q
ggsave(here::here("figures","PCAValues_mixedGenes.png"),plot = last_plot(),width = 7, height = 7)

