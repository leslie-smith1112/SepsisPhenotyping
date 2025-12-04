

## PLOT 

clus1all <- read.table(here::here("results","DECounted","Cluster1_DEALL_genes.tsv"), header = TRUE)
clus2all <- read.table(here::here("results","DECounted","Cluster2_DEALL_genes.tsv"), header = TRUE)
clus3all <- read.table(here::here("results","DECounted","Cluster3_DEALL_genes.tsv"), header = TRUE)
clus4all <- read.table(here::here("results","DECounted","Cluster4_DEALL_genes.tsv"), header = TRUE)
length(unique(c(clus1all$Symbol, clus2all$Symbol, clus3all$Symbol, clus4all$Symbol)))



## Datframe for colos
cluscolors <- data.frame("Cluster" = c("Cluster1Down", "Cluster1Up", "Cluster2Up","Cluster2Down","Cluster3Up","Cluster3Down","Cluster4Up", "Cluster4Down"),
                         "Colors" = c("#770059","#A80680" ,"#F55403", "#AF3D04" ,"#9381FF","#6054AB","#8AC926","#567E18"))

library(ggplot2)
clus1 <-  read.table(here::here("results","DifferentialExpression","Clus1VSAverage.tsv"), header = TRUE) |> dplyr::filter(Symbol %in% clus1all$Symbol)
head(clus1)
dim(clus1)
allDEGenes <- clus1
dim(allDEGenes)
#clus1 <- readclus1#clus1 <- readRDS(here::here("data","DifferentialExpression","Clus1VSAverage.rds"))
head(dat)
dat_or <- clus1[order(abs(clus1$logFC), decreasing = TRUE),]
head(dat_or)
dat_cut <- dat_or[1:50,]
dat_cut_or <- dat_cut[order(dat_cut$logFC, decreasing = TRUE),]
dat_cut_or$Symbol <- factor(dat_cut_or$Symbol, levels = dat_cut_or$Symbol)
dat_cut_or <- dat_cut_or |> dplyr::mutate("Direction" = ifelse(logFC <0, "Down", "Up" ))
dat_cut_or$Direction <- as.factor(dat_cut_or$Direction)
levels(dat_cut_or$Direction)
tail(dat_cut_or)
ggplot(dat_cut_or, aes(x = Symbol, y = logFC, fill = Direction, color = Direction),show.legend = FALSE) + 
  geom_bar(stat="identity", 
           show.legend = FALSE)+ 
  scale_fill_manual(values = c("#770059","#A80680")) +
  scale_color_manual(values = c("#770059","#A80680")) +
  theme_minimal() +
  theme(text = element_text(size = 23),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#+
scale_color_continuous(c(cluscolors$Color[cluscolors$Cluster == "Cluster1Down"], cluscolors$Color[cluscolors$Cluster == "Cluster1Up"]))


## PLOT 
clus2 <- read.table(here::here("results","DifferentialExpression","Clus2VSAverage.tsv"), header = TRUE) |> dplyr::filter(Symbol %in% clus2all$Symbol)
#clus2 <- readRDS(here::here("data","DifferentialExpression","Clus2VSAverage.rds"))
#dat <- tibble::rownames_to_column(clus2, "Genes")
dim(clus2)
allDEGenes <- rbind(allDEGenes, clus2)
dat_or <- clus2[order(abs(clus2$logFC), decreasing = TRUE),]
head(dat_or)
dat_cut <- dat_or[1:50,]
head(dat_cut)
range(dat_cut$logFC)
dat_cut_or <- dat_cut[order(dat_cut$logFC, decreasing = TRUE),]
dat_cut_or$Symbol <- factor(dat_cut_or$Symbol, levels = dat_cut_or$Symbol)
dat_cut_or <- dat_cut_or |> dplyr::mutate("Direction" = ifelse(logFC <0, "Down", "Up" ))
dat_cut_or$Direction <- as.factor(dat_cut_or$Direction)
ggplot(dat_cut_or, aes(x = Symbol, y = logFC,fill = Direction, color = Direction),show.legend = FALSE) + 
  geom_bar(stat="identity",  show.legend = FALSE)+ 
  scale_fill_manual(values = c("#AF3D04","#F55403")) +
  scale_color_manual(values = c("#AF3D04","#F55403")) +
  theme_minimal() +
  theme(text = element_text(size = 23),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


# Cluster3 Plot -----------------------------------------------------------
#clus3 <- readRDS(here::here("data","DifferentialExpression","Clus3VSAverage.rds"))
clus3 <- read.table(here::here("results","DifferentialExpression","Clus3VSAverage.tsv"), header = TRUE) |> dplyr::filter(Symbol %in% clus3all$Symbol)
dim(clus3)
allDEGenes <- rbind(allDEGenes, clus3)
dat_or <- clus3[order(abs(clus3$logFC), decreasing = TRUE),]
head(dat_or)
dat_cut <- dat_or[1:50,]
range(dat_cut$logFC)
dat_cut_or <- dat_cut[order(dat_cut$logFC, decreasing = TRUE),]
dat_cut_or$Symbol <- factor(dat_cut_or$Symbol, levels = dat_cut_or$Symbol)
dat_cut_or <- dat_cut_or |> dplyr::mutate("Direction" = ifelse(logFC <0, "Down", "Up" ))
dat_cut_or$Direction <- as.factor(dat_cut_or$Direction)
ggplot(dat_cut_or, aes(x = Symbol, y = logFC,fill = Direction, color = Direction),show.legend = FALSE) + 
  geom_bar(stat="identity",show.legend = FALSE)+ 
  scale_fill_manual(values = c("#6054AB","#9381FF")) +
  scale_color_manual(values = c("#6054AB","#9381FF")) +
  theme_minimal() +
  theme(text = element_text(size = 23),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


# Cluster 4 Plot ----------------------------------------------------------
clus4 <- read.table(here::here("results","DifferentialExpression","Clus4VSAverage.tsv"), header = TRUE) |> dplyr::filter(Symbol %in% clus4all$Symbol)
#clus4 <- readRDS(here::here("data","DifferentialExpression","Clus4VSAverage.rds"))
dim(clus4)
allDEGenes <- rbind(allDEGenes, clus4)
#write all DE genes to file:
write.table(allDEGenes, here::here("results","DifferentialExpression","AllDEGenes.tsv"), sep = '\t', col.names = TRUE, row.names = FALSE)
dat_or <- clus4[order(abs(clus4$logFC), decreasing = TRUE),]

head(dat_or)
dat_cut <- dat_or[1:50,]
range(dat_cut$logFC)
dat_cut_or <- dat_cut[order(dat_cut$logFC, decreasing = TRUE),]
dat_cut_or$Symbol <- factor(dat_cut_or$Symbol, levels = dat_cut_or$Symbol)
dat_cut_or <- dat_cut_or |> dplyr::mutate("Direction" = ifelse(logFC <0, "Down", "Up" ))
dat_cut_or$Direction <- as.factor(dat_cut_or$Direction)
ggplot(dat_cut_or, aes(x = Symbol, y = logFC,fill = Direction, color = Direction),show.legend = FALSE) + 
  geom_bar(stat="identity",  show.legend = FALSE)+ 
  scale_fill_manual(values = c("#567E18","#8AC926")) +
  scale_color_manual(values = c("#567E18","#8AC926")) +
  theme_minimal() +
  theme(text = element_text(size = 23),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


