
# differential gene  ------------------------------------------------------
clus1all <- read.table(here::here("results","differential_expression","Cluster1_DEALL_genes.tsv"), header = TRUE)
dim(clus1all)
length(unique(clus1all$Symbol[clus1all$adj.P.Val < 0.05]))
length(unique(clus1all$Symbol))
totab <- 
clus2all <- read.table(here::here("results","differential_expression","Cluster2_DEALL_genes.tsv"), header = TRUE)
length(unique(clus2all$Symbol))
length(unique(clus2all$Symbol[clus2all$adj.P.Val < 0.05]))
clus3all <- read.table(here::here("results","differential_expression","Cluster3_DEALL_genes.tsv"), header = TRUE)
length(unique(clus3all$Symbol[clus3all$adj.P.Val < 0.05]))
clus4all <- read.table(here::here("results","differential_expression","Cluster4_DEALL_genes.tsv"), header = TRUE)
length(unique(clus4all$Symbol))
length(unique(clus4all$Symbol[clus4all$adj.P.Val < 0.05]))
allgenes <- rbind(clus1all,clus2all,clus3all,clus4all)
clus2all |> 
  dplyr::distinct(Symbol, .keep_all = TRUE) |>  # Keep only the first row per SampleID
  dplyr::count(direction)



# read in all pathways ----------------------------------------------------

pathways <- read.table(here::here("results","KM","HumanBaseSignaturesKM","AllClusterPathwaysOrderedByClusterandModule.tsv"), header = TRUE)
pathways <- pathways[pathways$TERM_Q_VALUE <= 0.05,]
dim(pathways)
pathways <- dat # read on form pathwayscountsKM
pathcount <- table(pathways$CLUSTER[pathways$DIRECTION == "UP"], pathways$CLUSTER_NAME[pathways$DIRECTION == "UP"]) #17 down modules
apply(pathcount,1,sum)

# get counts of pathways 

df <- pathways |>
  dplyr::add_count(TERM_NAME, name = "All_Frequency")


df# add a summary string 
df2 <- df %>%
  # Step 1: For each TERM_NAME, if All_Frequency > 1 in any row, collect the CLUSTER info
  group_by(TERM_NAME) %>%
  mutate(
    Summary_Info = if (any(All_Frequency > 1)) {
      paste(unique(paste(CLUSTER, CLUSTER_NAME, DIRECTION, sep = ":")), collapse = "; ")
    } else {
      NA_character_
    }
  ) %>%
  dplyr::ungroup()
df2Order <- df2[order( df2$CLUSTER, df2$TERM_Q_VALUE ,df2$TERM_NAME),]
df2Order <- df2[order( df2$All_Frequency, df2$TERM_NAME, df2$CLUSTER, df2$CLUSTER_NAME ),]
head(df2Order)
#write.table(df2Order, here::here("final_files","AllPathwaysAddedSummaryColumnOrderedFrequency.tsv"), sep = "\t", col.names = TRUE, row.names = FALSE)

write.table(df2Order, here::here("results","KM","HumanBaseSignaturesKM","AllPathwaysAddedSummaryColumnOrderedFrequency.tsv"), sep = "\t", col.names = TRUE, row.names = FALSE)

mycut <- df |> dplyr::select(CLUSTER, TERM_NAME, DIRECTION)
head(mycut)

getoverlap <- function(clusterID){
  cat(" ## --------- DOING COMPARISON OF ", clusterID, " --------- ##\n" )
  cat("-------- Cluster 1 overlap --------\n")
  cat("DOWN")
  print(length(intersect(mycut$TERM_NAME[mycut$CLUSTER == clusterID & mycut$DIRECTION == "DOWN"], mycut$TERM_NAME[mycut$CLUSTER == "Cluster1" & mycut$DIRECTION == "DOWN"])))
  cat("UP")
  print(length(intersect(mycut$TERM_NAME[mycut$CLUSTER == clusterID & mycut$DIRECTION == "UP"], mycut$TERM_NAME[mycut$CLUSTER == "Cluster1" & mycut$DIRECTION == "UP"])))
  
  cat("-------- Cluster 2 overlap --------\n")
  cat("DOWN")
  print(length(intersect(mycut$TERM_NAME[mycut$CLUSTER == clusterID & mycut$DIRECTION == "DOWN"], mycut$TERM_NAME[mycut$CLUSTER == "Cluster2" & mycut$DIRECTION == "DOWN"])))
  cat("UP")
  print(length(intersect(mycut$TERM_NAME[mycut$CLUSTER == clusterID & mycut$DIRECTION == "UP"], mycut$TERM_NAME[mycut$CLUSTER == "Cluster2" & mycut$DIRECTION == "UP"])))
  
  cat("-------- Cluster 3 overlap --------\n")
  cat("DOWN")
  print(length(intersect(mycut$TERM_NAME[mycut$CLUSTER == clusterID & mycut$DIRECTION == "DOWN"], mycut$TERM_NAME[mycut$CLUSTER == "Cluster3" & mycut$DIRECTION == "DOWN"])))
  cat("UP")
  print(length(intersect(mycut$TERM_NAME[mycut$CLUSTER == clusterID & mycut$DIRECTION == "UP"], mycut$TERM_NAME[mycut$CLUSTER == "Cluster3" & mycut$DIRECTION == "UP"])))
  
  cat("-------- Cluster 4 overlap --------\n")
  cat("DOWN")
  print(length(intersect(mycut$TERM_NAME[mycut$CLUSTER == clusterID & mycut$DIRECTION == "DOWN"], mycut$TERM_NAME[mycut$CLUSTER == "Cluster4" & mycut$DIRECTION == "DOWN"])))
  cat("UP")
  print(length(intersect(mycut$TERM_NAME[mycut$CLUSTER == clusterID & mycut$DIRECTION == "UP"], mycut$TERM_NAME[mycut$CLUSTER == "Cluster4" & mycut$DIRECTION == "UP"])))
  
}
sink(here::here("results","KM","HumanBaseSignaturesKM","ClusterPathwayComparisons.log"), append = FALSE)
getoverlap("Cluster1")
getoverlap("Cluster2")
getoverlap("Cluster3")
getoverlap("Cluster4")
sink()

getoverlapopps <- function(clusterID){
  cat(" ## --------- DOING COMPARISON OF ", clusterID, " --------- ##\n" )
  cat("-------- Cluster 1 overlap --------\n")
  cat(paste0(clusterID, " DOWN", " Cluster 1 UP"))
  print(length(intersect(mycut$TERM_NAME[mycut$CLUSTER == clusterID & mycut$DIRECTION == "DOWN"], mycut$TERM_NAME[mycut$CLUSTER == "Cluster1" & mycut$DIRECTION == "UP"])))
  cat(paste0(clusterID, " UP", " Cluster 1 DOWN"))
  print(length(intersect(mycut$TERM_NAME[mycut$CLUSTER == clusterID & mycut$DIRECTION == "UP"], mycut$TERM_NAME[mycut$CLUSTER == "Cluster1" & mycut$DIRECTION == "DOWN"])))
  
  cat("-------- Cluster 2 overlap --------\n")
  cat(paste0(clusterID, " DOWN", " Cluster 2 UP"))
  print(length(intersect(mycut$TERM_NAME[mycut$CLUSTER == clusterID & mycut$DIRECTION == "DOWN"], mycut$TERM_NAME[mycut$CLUSTER == "Cluster2" & mycut$DIRECTION == "UP"])))
  cat(paste0(clusterID, " UP", " Cluster 2 DOWN"))
  print(length(intersect(mycut$TERM_NAME[mycut$CLUSTER == clusterID & mycut$DIRECTION == "UP"], mycut$TERM_NAME[mycut$CLUSTER == "Cluster2" & mycut$DIRECTION == "DOWN"])))
  
  cat("-------- Cluster 3 overlap --------\n")
  cat(paste0(clusterID, " DOWN", " Cluster 3 UP"))
  print(length(intersect(mycut$TERM_NAME[mycut$CLUSTER == clusterID & mycut$DIRECTION == "DOWN"], mycut$TERM_NAME[mycut$CLUSTER == "Cluster3" & mycut$DIRECTION == "UP"])))
  cat(paste0(clusterID, " UP", " Cluster 3 DOWN"))
  print(length(intersect(mycut$TERM_NAME[mycut$CLUSTER == clusterID & mycut$DIRECTION == "UP"], mycut$TERM_NAME[mycut$CLUSTER == "Cluster3" & mycut$DIRECTION == "DOWN"])))
  
  cat("-------- Cluster 4 overlap --------\n")
  cat(paste0(clusterID, " DOWN", " Cluster 4 UP"))
  print(length(intersect(mycut$TERM_NAME[mycut$CLUSTER == clusterID & mycut$DIRECTION == "DOWN"], mycut$TERM_NAME[mycut$CLUSTER == "Cluster4" & mycut$DIRECTION == "UP"])))
  cat(paste0(clusterID, " UP", " Cluster 4 DOWN"))
  print(length(intersect(mycut$TERM_NAME[mycut$CLUSTER == clusterID & mycut$DIRECTION == "UP"], mycut$TERM_NAME[mycut$CLUSTER == "Cluster4" & mycut$DIRECTION == "DOWN"])))
  
}



sink(here::here("results","KM","HumanBaseSignaturesKM","ClusterPathwayComparisonsOppositeDirections.log"), append = FALSE)
getoverlapopps("Cluster1")
getoverlapopps("Cluster2")
getoverlapopps("Cluster3")
getoverlapopps("Cluster4")
sink()




# Module specific overlaps:
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)

# Step 1: Create a named list of TERM_NAMEs by CLUSTER + CLUSTER_NAME
term_lists <- df2Order %>%
  group_by(GroupID = paste(CLUSTER, CLUSTER_NAME, DIRECTION, sep = "_")) |>
  summarise(terms = list(unique(TERM_NAME))) |>
  deframe()  # Converts to named list


group_combos <- expand.grid(Group1 = names(term_lists), Group2 = names(term_lists), stringsAsFactors = FALSE)

# Step 3: Calculate intersections
intersection_df <- group_combos %>%
  mutate(
    Shared_TERMS = map2_int(Group1, Group2, ~length(intersect(term_lists[[.x]], term_lists[[.y]]))),
    Intersecting_Terms = map2_chr(Group1, Group2, ~paste(intersect(term_lists[[.x]], term_lists[[.y]]), collapse = "; "))
  )



# Optional: Filter out self-comparisons or empty intersections
intersection_df_filtered <- intersection_df %>%
  filter(Group1 != Group2)

write.table(intersection_df_filtered, here::here("results","KM","HumanBaseSignaturesKM","ModuleSpecificClusterOverlapNONonSigs.tsv"), sep = '\t', col.names = TRUE, row.names = FALSE)
tt <- intersection_df_filtered |> filter(Shared_TERMS > 29)
tm <- tt[order(tt$Shared_TERMS),]
write.table(tm, here::here("results","KM","HumanBaseSignaturesKM","ModuleSpecificClusterOverlapTop46.tsv"), sep = '\t', col.names = TRUE, row.names = FALSE)


