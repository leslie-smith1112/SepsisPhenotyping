## Reference: https://davetang.github.io/muse/pheatmap.html#Heatmap_annotations
# Colors: https://colorhunt.co/

# Functions to save the pheatmap when running in terminal

save_pheatmap_pdf <- function(x, filename, width=12, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_png <- function(x, filename, width=12, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  png(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}


#temp
# 1.0 Read in data ----

# Read in expression matrix ----
expr <- readRDS(here::here("data","disease_expression.rds"))

# Read in genes used to define each cluster
# model_coeff <- readr::read_tsv(here::here("Multinomial_Model_Results","K5_Cluster_LASSO","Model_coefficients.tsv"))
# head(model_coeff)
# dim(model_coeff)

# model_coeff <- readr::read_tsv(here::here("results","Multinomial_Model_Results","PAM_K5","Gene_Signature.txt"))
# head(model_coeff)
# dim(model_coeff)
# 
# model_coeff <- readr::read_tsv(here::here("results","Multinomial_Model_Results","PAM_K5","Model_coeffifients_TopGenes.tsv"))
# model_coeff <- readr::read_tsv(here::here("results","Multinomial_Model_Results","PAM_K5","Model_coeffifients_TopGenes.tsv"))
dat <- readr::read_tsv(here::here("results","Multinomial_Model_Results","PAM_K4","Model_coefficients.tsv"))

dim(dat)
dat_non <- dat[dat$estimate != 0,]
length(unique(dat_non$term))
dat_non <- dat_non[dat_non$term != "(Intercept)",]
dat_sorted <- dat_non[order(dat_non$estimate,decreasing = TRUE),]
length(unique(dat_sorted$term))

# Adding in annotation of lipid genes -------------------------------------

gene_dat <- data.frame("Gene" = unique(dat_non$term), "Source" = rep("Variable",length(unique(dat_non$term))))

genes <- data.frame( Gene = c("CYP51A1","CYP4F2","LDLR","DHCR24","DHCR7","MSMO1","HMGCR","PCSK9","FDFT1","SQLE","ALOX5", "LSS","LCAT","LBR","CYP4Z1",
                              "KCNH7", "CYP4X1","CYP46A1","APOA1","EBP","CYP4V2", "CYP39A1","SCARB1","HSD17B7","CYP4F22","ALOX15","ABCG1","NSDHL", "CYP4F12","ALOX12",
                              "ABCA1","CYP4F11", "LOX","LIPE","TM7SF2","CYP4F8","PTGS1","LPA","LBR","CYP4F3","APOB","LBP","CETP","CYP4B1","CYP4A11",
                              "FDFT1","LCAT","CYP4A22","PLTP","PTGS2"), "Source" = "Lipid Genes")

gb_database <- readr::read_tsv(here::here("GO_Terms", "QuickGO-annotations.tsv"))

gb <- data.frame("Genes" = unique(gb_database$SYMBOL), "Source" = "GO DB")

gene_dat$Source[gene_dat$Gene %in% genes$Gene] <- "Lipid"
gene_dat$Source[gene_dat$Gene %in% gb$Genes] <- "GO"



expr_cut <- expr[rownames(expr) %in% gene_dat$Gene,]
dim(expr_cut)
# 2.0 Prep metadata for annotations ----

meta <- readRDS(here::here("data","Metadata_Filtered_withK.rds"))
meta$Gender[meta$Gender %in% c("H", "other", "Transgender")] <- "Other"
meta$Gender[!meta$Gender %in% c("Other", "Male", "Female")] <- "NoGenderReported"
table(meta$Gender)
meta_cut <- meta |> dplyr::select(`Sample ID`, Gender, AgeSummary, `Disease Simplified`, K4)
table(meta_cut$`Disease Simplified`)

# Redo levels
meta_cut$`Disease Simplified` <- factor(meta_cut$`Disease Simplified`, levels = c("Sepsis","Sepsis - Shock"))
head(meta_cut)
# Create Levels 
meta_cut$K4 <- as.factor(meta_cut$K4)
meta_cut$Gender <- as.factor(meta_cut$Gender)
meta_cut$AgeSummary <- as.factor(meta_cut$AgeSummary)
meta_cut$`Disease Simplified` <- as.factor(meta_cut$`Disease Simplified`)
head(meta_cut)
meta_cut <- tibble::column_to_rownames(meta_cut, "Sample ID")


# 3.0 Sort by Cluster → Disease → Age -> Gender ----
annotations <- meta_cut[order(meta_cut$K4, meta_cut$`Disease Simplified`, meta_cut$AgeSummary, meta_cut$Gender), ]
head(annotations)
colnames(annotations) <- c("Gender", "Age", "Disease", "Cluster ID")
annotations[1:5,1:4]
colnames(meta_cut) <- c("Gender", "Age", "Disease", "Cluster ID")
meta_cut <- meta_cut[,1:4]
head(meta_cut)
used[1:5,1:2]
row_annot <- gene_dat
row_annot$Source <- as.factor(row_annot$Source)
row_annot <- tibble::column_to_rownames(row_annot, "Gene")
head(row_annot)
# 3.1 Get columns and rows in the order we want
data_matrix <- expr_cut[, row.names(annotations)]
data_matrix[1:5,1:5]
dim(data_matrix)
ann_colors_col <- list(
  `Cluster ID` = c(Cluster1 = "#640D5F", Cluster2 = "#ffcfef", Cluster3 = "#DA498D", Cluster4 = "#A35C7A", Cluster5 = "#c5baff"),
  Gender = c(Female = "#4B5949", Male = "#B2C9AD", NoGenderReported = "#E5E3D4", Other = "#123458"),
  Age = c(Young = "#FADA7A", `Middle-Age` = "#81BFDA",Older = "#c4f4f9",NoAgeReported = "#E5E3D4"),
  Disease = c(Sepsis ="#f29F38", `Sepsis - Shock` = "#AF4459"),
  Source = c(GO = "#F564A9", Lipid = "#123458", Variable = "#CAE8BD")
)
data_matrix[1:5,1:5]

data_mat_cut <- data_matrix[rownames(data_matrix) %in% dat_sorted$term[1:20],]
dim(data_mat_cut)
data_mat_cut[1:5,1:5]

cluscolors <- data.frame("Cluster" = c("Cluster1", "Cluster2", "Cluster3", "Cluster4"), "Color" = c("#770059","#FF5400","#9381FF","#8AC926"))

range(expr_cut)
range(data_matrix)
# breakies = seq(-40, 100, 5)
# 
breakies <- seq(min(data_matrix), max(data_matrix), length.out = 100)
breakies <- seq(-5, 55, length.out = 75)  # Define color breaks from -10 to 100
colors <- colorRampPalette(c("blue", "white", "red"))(length(breakies) ) 
#colors <- colorRampPalette(c("blue", "white", "red"))(100)
# xx <- pheatmap::pheatmap(data_matrix, annotation_col = meta_cut, annotation_row = row_annotations, fontsize = 6, cluster_rows = FALSE, cluster_cols = FALSE,
#                          show_rownames = FALSE, show_colnames = FALSE, main = "Gene x Sample Correlation",
#                          annotation_colors = ann_colors_col, breaks = breakies,
#                          color = colors)
xx <- pheatmap::pheatmap(data_matrix, annotation_col = annotations, annotation_row = row_annot, fontsize = 5, cluster_rows = TRUE, cluster_cols = FALSE,
                         show_rownames = FALSE, show_colnames = FALSE, main = "Cluster Prediction Signature",
                         annotation_colors = ann_colors_col, breaks = breakies,
                         color = colors)
xx
save_pheatmap_pdf(xx, here::here("results","Multinomial_Model_Results","PAM_K4","GenexSampleLASSO_Heatmap2.pdf"))
save_pheatmap_png(xx, here::here("results","Multinomial_Model_Results","PAM_K4","GenexSampleLASSO_Heatmap2.png"))

xx <- pheatmap::pheatmap(data_mat_cut , annotation_col = annotations, fontsize = 7, cluster_rows = TRUE, cluster_cols = FALSE,
                         show_rownames = TRUE, show_colnames = FALSE, main = "Top Gene in Classification x Sample Correlation",
                         annotation_colors = ann_colors_col, breaks = breakies,
                         color = colors)
xx




dist_mat <- dist(t(data_mat_cut), method = 'euclidean')
hclust_avg <- hclust(dist_mat, method = 'average')
plot(hclust_avg)
cut_avg <- cutree(hclust_avg, k = 5)
suppressPackageStartupMessages(library(dendextend))
avg_dend_obj <- as.dendrogram(hclust_avg)
plot(hclust_avg)
rect.hclust(hclust_avg , k = 5, border = 2:6)
abline(h = 3, col = 'red')
avg_col_dend <- color_branches(avg_dend_obj, h = 5)
plot(avg_col_dend)




