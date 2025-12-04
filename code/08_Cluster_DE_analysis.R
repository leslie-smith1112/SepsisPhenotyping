## Clsuter analysis:

# 1.0 Read in all cluster DE for Cluster 1 :----

## listing where
list.files(here::here("results","differential_expression"))

readClus <- function(fileName){
  temp <- readRDS(here::here("results","differential_expression", fileName))
  return(temp)
}

de_summary <- function(dat1,dat2,dat3,clusterNum){
  
  all_clus <- rbind(dat1,dat2,dat3)
  cut_clus <- all_clus[all_clus$logFC != 0 ,] # get rid of genes with no logFC
  
  cut_clus$direction <- ifelse(cut_clus$logFC > 0, "up","down")
      
  # get value for whether a gene is differentially expressed the same way when compared to all other clusters
  gene_summary <- cut_clus |>
    dplyr::group_by(Symbol, direction) |>
    dplyr::summarize(count = dplyr::n(), .groups = "drop")
  table(gene_summary$count) #839 genes in here that are same in all clusters.
  message("Head:")
  print(head(gene_summary[gene_summary$count == 3,]))
  message("Tail")
  print(tail(gene_summary[gene_summary$count == 3,]))
  message("Table of counts:")
  print(table(gene_summary$count))
  to_write <- gene_summary[gene_summary$count == 3,]
  trimmed <- cut_clus[cut_clus$Symbol %in% to_write$Symbol,]
  #to_write <- merge(to_write, cut_clus, by = c("Symbol", "direction"))
  trim_or <- trimmed[order(abs(trimmed$logFC), decreasing = TRUE),]
  #print(to_write[1:30,])
  write.table(trim_or,  here::here("results","differential_expression", paste0(clusterNum,"_DEALL_genes.tsv")), sep = "\t", col.names = TRUE, row.names = FALSE)
  return(gene_summary)
}



clus12 <- readClus("Clust1VSClust2.tsv")  |> tibble::rownames_to_column("Symbol")
clus13 <- readClus("Clust1VSClust3.tsv")  |> tibble::rownames_to_column("Symbol")
clus14 <- readClus("Clust1VSClust4.tsv")  |> tibble::rownames_to_column("Symbol")


clus1_summ<- de_summary(clus12,clus13,clus14,"Cluster1") #839 # findign genes for a cluster that are up or down regulated in all cluster compares


# Cluster 2:--------
clus21 <- clus12
# we have to invert signs in this analysis because this was done with clsuter 1 as baseline
clus21$logFC <- -clus21$logFC 
head(clus21)
clus23 <- readClus("Clust2VSClust3.tsv") |> tibble::rownames_to_column("Symbol")
clus24 <- readClus("Clust2VSClust4.tsv") |> tibble::rownames_to_column("Symbol")

clus2_summ <- de_summary(clus21, clus23, clus24, "Cluster2") #606 genes


# Cluster 3:--------
# invert clus 1 v 3
clus31 <- clus13
clus31$logFC <- -clus31$logFC
# invert clus 2 v 3
clus32 <- clus23
clus32$logFC <- -clus32$logFC
clus34 <- readClus("Clust3VSClust4.tsv") |> tibble::rownames_to_column("Symbol")
clus3_summ <- de_summary(clus31, clus32, clus34, "Cluster3") #802 genes

# Cluster 4:--------
# invert 4 v 1 
clus41 <- clus14
clus41$logFC <- -clus41$logFC
#invert 4 v 2 
clus42 <- clus24
clus42$logFC <- -clus42$logFC
#invert 4 v 3
clus43 <- clus34
clus43$logFC <- -clus43$logFC
clus4_summ <- de_summary(clus41, clus42, clus43, "Cluster4") # 224


library(clusterProfiler)
library(org.Hs.eg.db)

# Initial investigatin of pathway hits with GSEA - not run for final analysis
# get GO 
myFUNC <- function(in_dat, file_name = NULL){ 
  #take the genes that are consistently up/down regualted in this cluster
  clus_dat <- in_dat[in_dat$count == 3,]
  print(table(clus_dat$count))
  # specify the up and down regulated genes
  up <- clus_dat[clus_dat$direction == "up",]
  down <- clus_dat[clus_dat$direction == "down",]
  dim(down)
  dim(up)
  up_proc <- enrichGO(gene =up$Symbol,
                              OrgDb = org.Hs.eg.db,
                              keyType = "SYMBOL",
                              ont = "ALL",
                              pAdjustMethod = "BH")
  write.table(up_proc, here::here("results","GSEA", paste0(file_name,"_UPREG.tsv")), sep = "\t", col.names = TRUE, row.names = FALSE)
  
  down_proc <- enrichGO(gene =down$Symbol,
                      OrgDb = org.Hs.eg.db,
                      keyType = "SYMBOL",
                      ont = "ALL",
                      pAdjustMethod = "BH")
  write.table(down_proc, here::here("results","GSEA", paste0(file_name,"_DOWNREG.tsv")), sep = "\t", col.names = TRUE, row.names = FALSE)
  
  clus_dat_proc <- enrichGO(gene =clus_dat$Symbol,
                        OrgDb = org.Hs.eg.db,
                        keyType = "SYMBOL",
                        ont = "ALL",
                        pAdjustMethod = "BH")
  write.table(clus_dat_proc, here::here("results","GSEA", paste0(file_name,"_ALLREG.tsv")), sep = "\t", col.names = TRUE, row.names = FALSE)
  
  
}

myFUNC(clus1_summ, "Cluster1")
myFUNC(clus2_summ, "Cluster2")
myFUNC(clus3_summ, "Cluster3")
myFUNC(clus4_summ, "Cluster4")





