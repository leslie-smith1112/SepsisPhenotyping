# Resources:
#https://erikaduan.github.io/posts/2021-01-02-volcano-plots-with-ggplot2/
#https://biostatsquid.com/volcano-plots-r-tutorial/


# Load libraries ----------------------------------------------------------
library(ggsurvfit)
library(survival)
library(survminer)
library(survRM2)
library(ggplot2)
library(svglite)

# Read in Metadata --------------------------------------------------------
meta <- readRDS(here::here("data","disease_metadata.rds"))
dim(meta)
colnames(meta)
meta_survive <- meta[meta$Survival %in% c("Survived", "Died") & !(is.na(meta$`Time to Event`)),]
dim(meta_survive)

# Survival Analysis -------------------------------------------------------
sur <- meta_survive |> dplyr::select(MolecularSubtype, TimetoEventEdited, DaySurvivalEdited)
head(sur)
sur$event <- 0
sur$event[sur$DaySurvivalEdited == "Died"] <- 1

# # Survival Object
surv_object <- Surv(time = sur$TimetoEventEdited ,event = sur$event)
fit <- survfit(surv_object ~ MolecularSubtype, data = sur)
wlrt(formula=Surv(event_time,event_status)~group,
     data=sim_data,
     method="mw",
     s_star = 0.5)


#  RMSt  -------------------------------------------------------
# 1 VS all RMST -----------------------------------------------------------

getRMST <- function(arm1, surc1){
  surc1$arm <- 0
  surc1$arm[surc1$MolecularSubtype != arm1] <- 1
  time=surc1$TimetoEventEdited
  status=surc1$event
  arm=as.numeric(surc1$arm)
  table(surc1$arm)
  sink(here::here("results","survival_analysis","RMST1VSAll.log"), append = TRUE)
  cat("---------------- ARM0: ",toupper(arm1)," VS OTHER ----------------")
  a=rmst2(time, status, arm)
  print(a)
  sink()
  plot(a, xlab="Years", ylab="Probability", density=60)
  
}
#C1
getRMST("Cluster1",sur)
p.adjust(0.49, method = "bonferroni", n = 4)
#C2
getRMST("Cluster2",sur)
p.adjust(0.067, method = "bonferroni", n = 4)
#C3
getRMST("Cluster3",sur)
p.adjust(0.067, method = "bonferroni", n = 4)
#C4
getRMST("Cluster4",sur)
p.adjust(0.009, method = "bonferroni", n = 4)


# Survival Plot -----------------------------------------------------------
ggsurvplot(fit,
           data = sur,
           pval = TRUE,             # Show p-value of log-rank test
           conf.int = FALSE,  
           xlim = c(0, 28),# Show confidence intervals
           risk.table = TRUE,# Show number at risk table
           break.x.by = 4, 
           surv.median.line = "hv", # Add median survival lines
           palette = c("#770059", "#F55403", "#9381FF", "#8AC926" ),       # Color palette
           legend.title = "Cluster",
           legend.labs = levels(sur$MolecularSubtype)) 




# Violin Plot for Age  ----------------------------------------------------
library(ggplot2)
meta_age <- meta[!(is.na(meta$Age)),]
head(meta_age)
dim(meta_age)
table(meta_age$Age)
age_cut <- meta_age |> dplyr::select(MolecularSubtype, Age, Dataset)
again <- age_cut[!(age_cut$Age %in% c("20-39","40-59","60-79", "> 80")) & !(is.na(age_cut$Age)),]
dim(again)
again$Age<- as.numeric(again$Age)

ggplot(again,aes(x = MolecularSubtype, y = Age, fill = MolecularSubtype)) +
  geom_violin(trim = FALSE) +
  scale_fill_manual(values = c("#770059", "#F55403", "#9381FF", "#8AC926" )) +
  geom_boxplot(width = 0.1, outlier.shape = NA, color = "black") +
  theme_classic() + geom_hline(yintercept = median(again$Age, na.rm = TRUE), 
                               linetype = "dashed", color = "black", linewidth = 0.5)

stat_summary(fun = median, geom = "point", 
             shape = 23, size = 3, fill = "white", color = "black") +
  geom_jitter(width = 0.2, size = 1, alpha = 0.4)




# Volcano Plot ------------------------------------------------------------

# Differential Expression
library(limma)
#library(edgeR)

# 1.0 Read in data ----
expr_dat <- readRDS(here::here("data","disease_expression.rds"))

#meta<- readRDS(here::here("data","Metadata_Filtered_withK.rds")) # made in the differnetial expression script
survive <- meta[meta$Survival %in% c("Died","Survived"),]
dim(survive) #(n = 1710)

# 2.0 Make sure the cluster sample IDs are in the same order as the colnames in expr_dat ----
expr_cut <- expr_dat[,colnames(expr_dat) %in% survive$`Sample ID`]
dim(expr_cut)
meta <- survive[match(colnames(expr_cut), survive$`Sample ID`), ]
dim(meta)

all(colnames(expr_cut) == meta$`Sample ID`)
dim(meta)

# 3.0 Create design matrix for DEs for K5

# ---- with K = 4 ----- # 
meta$Survival <- as.factor(meta$Survival)
levels(meta$Survival) # levels: [1] "Died"     "Survived"
design <- model.matrix(~Survival, data = meta ) #direct comparison of clusters - do not want intercept
# colnames(design) <- "Survival"
# head(design)



# 4.0 Fit the linear model ----

fit <- lmFit(expr_cut, design)
fit <- eBayes(fit, trend = TRUE, robust = TRUE)

# 4.1 Create contrast matrices for each cluster ----
## 4.1.2 Cluster 1 Comparisons ----
# Average Across all
# contrast.matrix <- makeContrasts(Survived - Died, levels=design)
# fit2 <- contrasts.fit(fit, contrast.matrix)
# 
# fit2 <- eBayes(fit2, trend = TRUE, robust = TRUE)


#res <- topTable(fit2, adjust.method = "BH", number = Inf)
res <- topTable(fit, adjust.method = "BH", number = Inf)
head(res)
res_edit  <- res |> mutate(fold_change = 2^logFC)
res_edit <- res_edit |> mutate(gene_type = dplyr::case_when(log2(fold_change) > 0 & adj.P.Val <= 0.05 ~ "up",
                                                            log2(fold_change) <0 & adj.P.Val <= 0.05 ~ "down",
                                                            TRUE ~ "ns"))

head(res_edit)
table(res_edit$gene_type)
res_edit <- res_edit[order(res_edit$gene_type, res_edit$logFC, decreasing = TRUE),]

saveRDS(res_edit, here::here("results", "survival_analysis","SurvivalDE.rds"))
res_edit <- tibble::rownames_to_column(res_edit, "Symbol")
head(res_edit)

write.table(res_edit, here::here("results", "survival_analysis","Survival_DE.tsv"), sep = '\t', col.names = TRUE, row.names = FALSE)


# Initially we had up/down in reference to survivors, group requests direction reference patients who died
res_edit$gene_type_edited <- NA
res_edit$gene_type_edited[res_edit$gene_type == "up"] <- "down"
res_edit$gene_type_edited[res_edit$gene_type == "down"] <- "up"
res_edit$gene_type_edited[res_edit$gene_type == "ns"] <- "ns"
temp <- res_edit[order(abs(log2(res_edit$fold_change)),decreasing = TRUE),]
who <- res_edit[order(abs(res_edit$logFC),decreasing = TRUE),]

# Genes to label in plot
upgenes <- temp[temp$Symbol %in% c("IL1R2",   "LCN2",  "CD177",  "MMP9","OLFM4","S100A12", "LTF","ARG1", "CA1", 
                                   "S100A8", "S100P",   "STOM",  "TTK", "TCN1","SAMSN1"," AREG","HMGB2", 
                                   "PPBP","CEACAM8","RETN","VNN1","ALOX5AP","PFKFB2", "CDKN3","RANBP9","GCA", "MYL6"),]
#Down in Died
downgenes <- temp[temp$Symbol %in% c( "CTSS","CD74", "ITM2B","HLA-DRA","FGD3","HLA-E","LITAF","FGL2","FFAR2","TGFBI",
                                      "GIMAP4", "CX3CR1","NAGK", "MX1", "HLA-DPA1", "HLA-C", "TMSB10","RGS2","GBP5",
                                      "GIMAP4"),]
# Sqwapping the values so that we have the reference to Died as faheem asked (I.E up is up in died )
head(res_edit)

allsig <- rbind(upgenes,downgenes)
allsig
cols <- c("up" = "red", "down" = "blue", "ns" = "grey") 
colscol <- c("up" = "red", "down" = "blue", "ns" = "#6b6a6a") 
sizes <- c("up" = 2, "down" = 2, "ns" = 1) 
alphas <- c("up" = 0.5, "down" = 0.5, "ns" = 0.5)
library(ggrepel)
ggplot(res_edit,aes(x = -log2(fold_change),
                    y = -log10(adj.P.Val),
                    fill = gene_type_edited,    
                    size = gene_type_edited,
                    alpha = gene_type_edited,
                    color = gene_type_edited)) + 
  geom_point(shape = 21 # Specify shape and colour as fixed local parameters    
  ) +
  geom_point(data = upgenes,
             shape = 21,
             size = 2, 
             fill = "#760202",
             colour = "#760202") + 
  geom_point(data = downgenes,
             shape = 21,
             size = 2, 
             fill = "navy",
             colour = "navy") + 
  geom_text_repel(data = allsig, # Add labels last to appear as the top layer  
                  aes(label = Symbol), size = 6, color = "black")+
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") + 
  geom_vline(xintercept = c(log2(0.5), log2(2)),
             linetype = "dashed") +
  scale_x_continuous(breaks = c(seq(-17, 17, 2)), # Modify x-axis tick intervals    
                     limits = c(-17, 17)) +
  scale_fill_manual(values = cols) +
  #scale_colour_manual(values = cols) +# Modify point colour
  scale_size_manual(values = sizes) + # Modify point size
  scale_alpha_manual(values = alphas) +
  scale_colour_manual(values = colscol) + theme_minimal()  + 
  theme(axis.text=element_text(size=20)) +
  labs(title = "Gene expression changes in Survival",
       x = "log2(fold change)",
       y = "-log10(adjusted P-value)")

