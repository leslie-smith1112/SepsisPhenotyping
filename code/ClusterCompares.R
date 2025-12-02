source(here::here("code","SilouetteScore.R"))


# Function Definition for Clustering -------------------------------------

Cluster <- function(cluster_type, norm_matrix, title){
  message("STARTING", toupper(cluster_type),"CLUSTERING","\n")
  if(cluster_type == "km"){
    results = ConsensusClusterPlus(norm_matrix ,maxK=7,reps = 1000,
                                 title=title,clusterAlg=cluster_type, distance="euclidean",seed=1262118388.71279,plot="png")
  }else{
    results = ConsensusClusterPlus(norm_matrix ,maxK=7,reps = 1000,
                                 title=title,clusterAlg="hc", distance="pearson",seed=1262118388.71279,plot="png")
  }
  # results = ConsensusClusterPlus(d,maxK=6,reps=50,pItem=0.8,pFeature=1,
  #                                + title=title,clusterAlg="hc",distance="pearson",seed=1262118388.71279,plot="png")
  dir.create(here::here("results","ClusterCompares", title))
  
  #get item consensus for the clusters
  icl = calcICL(results,title=title,plot="png")
  icl[["clusterConsensus"]]
  write.table(icl[["clusterConsensus"]],here::here("results","ClusterCompares",title,"ClusterConsensus_Assignments.tsv"), sep = "\t", col.names = TRUE, row.names = FALSE)
  
  #For every cluster tested, write the sample cluster assignments
  for(i in 2:7){#TODO is this a list?? change to lapply
    cluster_temp <- data.frame("Samples" = names(results[[i]][["consensusClass"]]), "ClusterID" = results[[i]][["consensusClass"]])
    cluster_temp$ClusterID <- paste0("Cluster", cluster_temp$ClusterID)
    write.table(cluster_temp,here::here("results","ClusterCompares",title,paste0("Sample_Cluster_Assignments_",i,".tsv")), sep = "\t", col.names = TRUE, row.names = FALSE)
    
  }
  message("FINISHED", toupper(cluster_type),"CLUSTERING","\n")
  message("STARTING SILOUETTE CALCULATIONS FOR", toupper(cluster_type),"CLUSTERING","\n")
  silouette_scores(swept_away, title, "normal")
  message("COMPLETED SILOUETTE CALCULATIONS FOR", toupper(cluster_type),"CLUSTERING","\n")
}



# 1.0 Read in data ----
expr_dat <- readRDS(here::here("data","expr_mixedGenes.rds"))
message("READ IN DATA\n")

# 2.0 normal clustering ----
library(ConsensusClusterPlus)

#row-wise median centering
swept_away <- sweep(expr_dat,1, apply(expr_dat,1,median,na.rm=T))
message("NORMALIZED DATA\n")
# Cluster with PAM, K Means and HC for comparisons

Cluster("pam",swept_away, "PAM_VALIDATE")
Cluster("hc",swept_away, "HC_VALIDATE")
Cluster("km",swept_away, "KM_VALIDATE")
  
 
#results = ConsensusClusterPlus(swept_away ,maxK=7,reps = 1000,
#                               title="HC",clusterAlg="hc", distance="pearson",seed=1262118388.71279,plot="png")
  
#dir.create(here::here("results","ClusterCompares", title))

#get item consensus for the clusters
#icl = calcICL(results,title=title,plot="png")
#icl[["clusterConsensus"]]
#write.table(icl[["clusterConsensus"]],here::here("results","ClusterCompares",title,"ClusterConsensus_Assignments.tsv"), sep = "\t", col.names = TRUE, row.names = FALSE)

#For every cluster tested, write the sample cluster assignments
#for(i in 2:7){#TODO is this a list?? change to lapply
#  cluster_temp <- data.frame("Samples" = names(results[[i]][["consensusClass"]]), "ClusterID" = results[[i]][["consensusClass"]])
#  cluster_temp$ClusterID <- paste0("Cluster", cluster_temp$ClusterID)
#  write.table(cluster_temp,here::here("results","ClusterCompares",title,paste0("Sample_Cluster_Assignments_",i,".tsv")), sep = "\t", col.names = TRUE, row.names = FALSE)
  
#}
#message("FINISHED", toupper(cluster_type),"CLUSTERING","\n")
#message("STARTING SILOUETTE CALCULATIONS FOR", toupper(cluster_type),"CLUSTERING","\n")
#silouette_scores(swept_away, title, "normal")
