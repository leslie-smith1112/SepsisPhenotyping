# Differential Expression
library(limma)
library(ggplot2)

# 1.0 Read in data ----
expr_dat <- readRDS(here::here("data","disease_expression.rds"))
meta <- readRDS(here::here("data","disease_metadata.rds"))


# 2.0 Make sure the cluster sample IDs are in the same order as the colnames in expr_dat ----
meta <- meta[match(colnames(expr_dat), meta$`Sample ID`), ]
all(colnames(expr_dat) == meta$`Sample ID`)
dim(meta)

# 3.0 Create design matrix for DEs for K5

# ---- with K = 4 ----- # 
meta$MolecularSubtype <- as.factor(meta$MolecularSubtype)
design <- model.matrix(~0 + MolecularSubtype, meta ) #direct comparison of clusters - do not want intercept
colnames(design) <- levels(meta$MolecularSubtype)
head(design)


# 4.0 Fit the linear model ----

fit <- lmFit(expr_dat, design)


# 4.1 Create contrast matrices for each cluster ----
## 4.1.2 Cluster 1 Comparisons ----
# Average Across all
contrast.matrix <- makeContrasts(Cluster1 - (Cluster2 + Cluster3 + Cluster4)/3, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)

fit2 <- eBayes(fit2, trend = TRUE, robust = TRUE)

head(fit2)
#res <- topTable(fit2, adjust.method = "BH", number = Inf)
res <- topTable(fit2, coef = 1, adjust.method = "BH", number = Inf)
head(res)

tail(res)
dim(res)
de <- res[res$adj.P.Val < 0.05 & abs(res$logFC) > 1, ]
de$Cluster <- "Cluster1VSAvgAll"
dim(de)
de <- tibble::rownames_to_column(de, "Symbol")
head(de)

write.table(de, here::here("results","differential_expression","Clus1VSAverage.tsv"), sep = '\t', col.names = TRUE, row.names = FALSE)


# dim(de)
# genes <- data.frame( Gene = c("CYP51A1","CYP4F2","LDLR","DHCR24","DHCR7","MSMO1","HMGCR","PCSK9","FDFT1","SQLE","ALOX5", "LSS","LCAT","LBR","CYP4Z1",
#                               "KCNH7", "CYP4X1","CYP46A1","APOA1","EBP","CYP4V2", "CYP39A1","SCARB1","HSD17B7","CYP4F22","ALOX15","ABCG1","NSDHL", "CYP4F12","ALOX12",
#                               "ABCA1","CYP4F11", "LOX","LIPE","TM7SF2","CYP4F8","PTGS1","LPA","LBR","CYP4F3","APOB","LBP","CETP","CYP4B1","CYP4A11",
#                               "FDFT1","LCAT","CYP4A22","PLTP","PTGS2"), "Source" = "Lipid Genes")
# dat_cut[dat_cut$Genes %in% 
#           genes$Gene,]


# First we will compare cluster 1 with all other clusters, then we will do vs specific clusters
#This should be any gene that is differnetiually epxressed between Cluster 1 and another clsuter
contrast.matrix <- makeContrasts(Cluster1 - Cluster2, Cluster1 - Cluster3, Cluster1 - Cluster4, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2, trend = TRUE, robust = TRUE)
# All variances
res <- topTable(fit2, coef=c(1,2,3), adjust.method = "BH", number=Inf)
head(res)
tail(res)
colnames(res) <- c("logFC_Clus1v2","logFC_Clus1v3","logFC_Clus1v4", "AveExpr","F", "P.Value", "adj.P.Val")
head(res)
de <- res[res$adj.P.Val < 0.05, ]
dim(de)
de$Cluster <- "Cluster1VSAllOthers"
dim(de)
de <- tibble::rownames_to_column(de, "Symbol")
head(de)
write.table(de, here::here("results","differential_expression","Clus1VSAllOthers.tsv"), sep = '\t', col.names = TRUE, row.names = FALSE)

#Cluster 1 VS 2 ----
res <- topTable(fit2, coef = 1, adjust.method = "BH", number = Inf)
head(res)
tail(res)
dim(res)
de <- res[res$adj.P.Val < 0.05 & abs(res$logFC) > 1, ] # can now filter by logFC too
dim(de)
head(de)
de$Cluster <- "Cluster1VSCluster2"
dim(de)
de <- tibble::rownames_to_column(de, "Symbol")
head(de)
write.table(de, here::here("results","differential_expression","Clus1VSClus2.tsv"), sep = '\t', col.names = TRUE, row.names = FALSE)
IM HERE

#Cluster 1 VS 3 ----
res <- topTable(fit2, coef=2, adjust.method = "BH", number=Inf)
head(res)
tail(res)
dim(res)
de <- res[res$adj.P.Val < 0.05 & abs(res$logFC) > 1, ] # can now filter by logFC too
dim(de)
head(de)
de$Cluster <- "Cluster1VSCluster3"
dim(de)
de <- tibble::rownames_to_column(de, "Symbol")
head(de)
write.table(de, here::here("results","differential_expression","Clus1VSClus3.tsv"), sep = '\t', col.names = TRUE, row.names = FALSE)

#Cluster 1 VS 4 ----
res <- topTable(fit2, coef=3,adjust.method = "BH", number=Inf)
head(res)
tail(res)
dim(res)
de <- res[res$adj.P.Val < 0.05 & abs(res$logFC) > 1, ] # can now filter by logFC too
dim(de)
head(de)
de$Cluster <- "Cluster1VSCluster4"
dim(de)
de <- tibble::rownames_to_column(de, "Symbol")
head(de)
write.table(de, here::here("results","differential_expression","Clus1VSClus4.tsv"), sep = '\t', col.names = TRUE, row.names = FALSE)



# 4.1.2 Cluster 2 -------------------------------------------------------------

# Average Across all
contrast.matrix <- makeContrasts(Cluster2 - (Cluster1 + Cluster3 + Cluster4)/3, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2, trend = TRUE, robust = TRUE)
head(fit2)
res <- topTable(fit2, coef = 1, adjust.method = "BH", number = Inf)
head(res)
tail(res)
de <- res[res$adj.P.Val < 0.05 & abs(res$logFC) > 1, ]
dim(de)
de$Cluster <- "Cluster2VSAvgAll"
dim(de)
de <- tibble::rownames_to_column(de, "Symbol")
head(de)
write.table(de, here::here("results","differential_expression","Clus2VSAverage.tsv"), sep = '\t', col.names = TRUE, row.names = FALSE)


contrast.matrix <- makeContrasts(Cluster2 - Cluster1, Cluster2 - Cluster3, Cluster2 - Cluster4, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2, trend = TRUE, robust = TRUE)
# All variances
res <- topTable(fit2, coef=c(1,2,3), adjust.method = "BH", number=Inf)
head(res)
tail(res)
colnames(res) <- c("logFC_Clus2v1","logFC_Clus2v3","logFC_Clus2v4", "AveExpr","F", "P.Value", "adj.P.Val")
head(res)
de <- res[res$adj.P.Val < 0.05, ]
dim(de)
de$Cluster <- "Cluster2VSAllOthers"
de <- tibble::rownames_to_column(de, "Symbol")
head(de)
write.table(de, here::here("results","differential_expression","Clus2VSOthers.tsv"), sep = '\t', col.names = TRUE, row.names = FALSE)


#Cluster 2 VS 3 ----
res <- topTable(fit2, coef=2, adjust.method = "BH", number=Inf)
head(res)
tail(res)
dim(res)
de <- res[res$adj.P.Val < 0.05 & abs(res$logFC) > 1, ] # can now filter by logFC too
dim(de)
head(de)
de$Cluster <- "Cluster2VSCluster3"
dim(de)
de <- tibble::rownames_to_column(de, "Symbol")
head(de)
write.table(de, here::here("results","differential_expression","Clus2VSClus3.tsv"), sep = '\t', col.names = TRUE, row.names = FALSE)


#Cluster 2 VS 4 ----
res <- topTable(fit2, coef=3, adjust.method = "BH", number=Inf)
head(res)
tail(res)
dim(res)
de <- res[res$adj.P.Val < 0.05 & abs(res$logFC) > 1, ] # can now filter by logFC too
dim(de)
de$Cluster <- "Cluster2VSCluster4"
dim(de)
de <- tibble::rownames_to_column(de, "Symbol")
head(de)
write.table(de, here::here("results","differential_expression","Clus2VSClus4.tsv"), sep = '\t', col.names = TRUE, row.names = FALSE)

# 4.1.3 Cluster 3 -------------------------------------------------------------

# Average Across all
contrast.matrix <- makeContrasts(Cluster3 - (Cluster1 + Cluster2 + Cluster4) /3, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2, trend = TRUE, robust = TRUE)
head(fit2)
res <- topTable(fit2, coef = 1, adjust.method = "BH", number = Inf)
head(res)
tail(res)
de <- res[res$adj.P.Val < 0.05 & abs(res$logFC) > 1, ]
dim(de)
de$Cluster <- "Cluster3VSAvgAll"
dim(de)
de <- tibble::rownames_to_column(de, "Symbol")
head(de)
write.table(de, here::here("results","differential_expression","Clus3VSAverage.tsv"), sep = '\t', col.names = TRUE, row.names = FALSE)


contrast.matrix <- makeContrasts(Cluster3 - Cluster1, Cluster3 - Cluster2, Cluster3 - Cluster4, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2, trend = TRUE, robust = TRUE)
# All variances
res <- topTable(fit2, coef=c(1,2,3), adjust.method = "BH", number=Inf)
head(res)
tail(res)
colnames(res) <- c("logFC_Clus3v1","logFC_Clus3v2","logFC_Clus3v4", "AveExpr","F", "P.Value", "adj.P.Val")
head(res)
de <- res[res$adj.P.Val < 0.05, ]
dim(de)
de$Cluster <- "Cluster3VSAllOthers"
dim(de)
de <- tibble::rownames_to_column(de, "Symbol")
head(de)
write.table(de, here::here("results","differential_expression","Clus3VSOthers.tsv"), sep = '\t', col.names = TRUE, row.names = FALSE)


#Cluster 3 VS 4 ----
res <- topTable(fit2, coef=3, adjust.method = "BH", number=Inf)
head(res)
tail(res)
dim(res)
de <- res[res$adj.P.Val < 0.05 & abs(res$logFC) > 1, ] # can now filter by logFC too
dim(de)
de$Cluster <- "Cluster3VSCluster4"
dim(de)
de <- tibble::rownames_to_column(de, "Symbol")
head(de)
write.table(de, here::here("results","differential_expression","Clus3VSClus4.tsv"), sep = '\t', col.names = TRUE, row.names = FALSE)


# 4.1.4 Cluster 4 -------------------------------------------------------------

# Average Across all
contrast.matrix <- makeContrasts(Cluster4 - (Cluster1 + Cluster2 + Cluster3)/3, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2, trend = TRUE, robust = TRUE)
head(fit2)
res <- topTable(fit2, coef = 1, adjust.method = "BH", number = Inf)
head(res)
tail(res)
de <- res[res$adj.P.Val < 0.05 & abs(res$logFC) > 1, ]
dim(de)
de$Cluster <- "Cluster4VSAvgAll"
dim(de)
de <- tibble::rownames_to_column(de, "Symbol")
head(de)
write.table(de, here::here("results","differential_expression","Clus4VSAverage.tsv"), sep = '\t', col.names = TRUE, row.names = FALSE)


contrast.matrix <- makeContrasts(Cluster4 - Cluster1, Cluster4 - Cluster2, Cluster4 - Cluster3, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2, trend = TRUE, robust = TRUE)
# All variances
res <- topTable(fit2, coef=c(1,2,3),adjust.method = "BH", number=Inf)
head(res)
tail(res)
colnames(res) <- c("logFC_Clus4v1","logFC_Clus4v2","logFC_Clus4v3", "AveExpr","F", "P.Value", "adj.P.Val")
head(res)
de <- res[res$adj.P.Val < 0.05, ]
dim(de)
de$Cluster <- "Cluster4VSAllOthers"
dim(de)
de <- tibble::rownames_to_column(de, "Symbol")
head(de)
write.table(de, here::here("results","differential_expression","Clus4VSOthers.tsv"), sep = '\t', col.names = TRUE, row.names = FALSE)

