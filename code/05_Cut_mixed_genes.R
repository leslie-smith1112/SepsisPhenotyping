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

# 2.0 Look at variable genes in the dataset'
gene_variability <- apply(expr_dat, 1, mad)
ranked_variances <- sort(gene_variability, decreasing = TRUE)
head(ranked_variances)
top_genes <- names(ranked_variances)[1:5000]
head(top_genes)



# 3.0 immune related genes
gb_database <- readr::read_tsv(here::here("GO_Terms", "QuickGO-annotations.tsv"))
gb_database[1:5,1:5]
length(unique(gb_database$SYMBOL))
length(intersect(rownames(expr_dat), gb_database$SYMBOL)) #1207 intersecting genes
expr_dat[1:5,1:5]

# 4.0 Lipid genes 
genes <- c("CYP51A1","CYP4F2","LDLR","DHCR24","DHCR7","MSMO1","HMGCR","PCSK9","FDFT1","SQLE","ALOX5", "LSS","LCAT","LBR","CYP4Z1",
           "KCNH7", "CYP4X1","CYP46A1","APOA1","EBP","CYP4V2", "CYP39A1","SCARB1","HSD17B7","CYP4F22","ALOX15","ABCG1","NSDHL", "CYP4F12","ALOX12",
           "ABCA1","CYP4F11", "LOX","LIPE","TM7SF2","CYP4F8","PTGS1","LPA","LBR","CYP4F3","APOB","LBP","CETP","CYP4B1","CYP4A11",
           "FDFT1","LCAT","CYP4A22","PLTP","PTGS2")


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
saveRDS(expr_mixed, here::here("data","disease_expression.rds"))as

# 6.0 Plot PCA
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
saveRDS(dtp_mixed, here::here("data","PCAValues_disease_expression.rds"))


# -- Plotting --
q <- ggplot(data=dtp_mixed, aes(x=PC1, y=PC2, color=Disease, shape = Technology)) + 
  geom_point(aes(alpha = 0.7)) + 
  scale_color_manual(values = c("#A59DA2","#a50a6c","#192dfa","#07d81f")) +
  theme_bw() +  
  ggtitle(paste0("GO + Lipid + Variance Genes with ",nrow(dtp_mixed)," Samples")) + 
  labs(y= "PC2", x = "PC1")
q
ggsave(here::here("figures","PCAValues_disease_expression.png"),plot = last_plot(),width = 7, height = 7)



