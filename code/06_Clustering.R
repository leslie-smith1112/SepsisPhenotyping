library(logr)
library(cluster)
library(stats)
library(png)
library(grid)
library(gridExtra)
require("cluster")

# Function Definition for Clustering -------------------------------------
silouette_scores <- function(normalized_mat, title, clus_type){
  # 1.3 Set matrix for cluster ----
  if(clus_type == "normal"){
    message("Creating distance matrix with swept gene-sample matrix")
    # for regular clustering:
    #swept_away <- sweep(expression,1, apply(expression,1,median,na.rm=T))
    #silhouette needs distance matrix - dist() goes by rows.
    dist_matrix <- dist(t(normalized_mat))
  }else if(clus_type == "correlation"){
    message("Creating distance matrix using sample-sample correlation")
    # for the correlation matrix: 
    correlated <- cor(expression, method = "spearman")
    #silhouette needs distance matrix
    dist_matrix <- as.dist(1 - correlated)
    
  }else{
    messasge("Please enter valid cluster type.")
    stop()
  }
  
  
  # 2.0 Silhouette scoring and plots ----
  
  
  # Read in cluster assignments for given title
  sink(here::here("results","cluster_compares",title,paste0("Sil_Scores.log")), append = TRUE)
  for(i in 2:7){
    assigns <- readr::read_tsv(here::here("results","cluster_compares",title,paste0("Sample_Cluster_Assignments_",i,".tsv")))
    #these are in the same order as the columns in expression - as they must be
    clusters <- as.numeric(substr(assigns$ClusterID, nchar(assigns$ClusterID),nchar(assigns$ClusterID)))
    k = max(clusters)
    clusters <- as.vector(clusters)
    colnames(assigns) <- c("Samples", paste0("K",i))
    
    
    ## 2.1 Create master matrix holding all K cluster assignments for this title ----
    
    if(i==2){
      master_dat <- assigns
    }else if(i>2){
      # this merge reorders the samples so we need to make sure the samples are in correct order if we use this
      master_dat <- merge(master_dat, assigns, by = "Samples") #should end up with a matrix that is 2251 x 4
    }
    
    
    ## 2.2  Calculate silhouette scores ----
    
    sil <- silhouette(clusters, dist_matrix)
    cat(paste0("---------- Cluster Num: ",i, " ---------- \n"))
    out_score <- summary(sil)
    print(summary(sil))
    cat("\n")
    
    
    # get silhouette plots
    plot1 <- factoextra::fviz_silhouette(sil)
    cp <- ggpubr::ggpar(plot1,subtitle = paste0("Clustering with Genes for K",i))
    c2 <- cp + ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank())# + ggpubr::grids(linetype = "dashed")
    
    
    ggplot2::ggsave(here::here("results","cluster_compares","plots",title, paste0("SilPlot_Genes_K",i,".png")), plot = c2, width = 8, height = 6, create.dir = TRUE)
  }
  sink()
  
  # * -- savve RDS -- * 
  saveRDS(master_dat, here::here("results","cluster_compares", title, "ClusterIDs.rds"))
  
  # 2.3 Master silhouette plot ----
  filenames <- here::here("results","cluster_compares","plots",title, paste0("SilPlot_Genes_K",2:7,".png"))
  
  piclist<-list()
  for(j in 1:6){piclist[[j]]<-readPNG(filenames[j])}
  
  n <- length(piclist)
  nCol <- floor(sqrt(n))
  
  holding <- arrangeGrob(rasterGrob(piclist[[1]]),rasterGrob(piclist[[2]]),
                         rasterGrob(piclist[[3]]),rasterGrob(piclist[[4]]),
                         rasterGrob(piclist[[5]]),rasterGrob(piclist[[6]]),
                         ncol=nCol)
  holding
  ggplot2::ggsave(here::here("results","cluster_compares","plots",title,"AllSilPlot_GeneKs.png"),holding,width=15,height=11, create.dir = TRUE)
}


Cluster <- function(cluster_type, norm_matrix, title){
  message("STARTING", toupper(cluster_type),"CLUSTERING","\n")
  if(cluster_type == "km"){
    results = ConsensusClusterPlus(norm_matrix ,maxK=7,reps = 1000, pItem=0.8,pFeature=1,
                                   title=title,clusterAlg=cluster_type, distance="euclidean",seed=1262118388.71279,plot="png")
  }else{
    results = ConsensusClusterPlus(norm_matrix ,maxK=7,reps = 1000, pItem=0.8,pFeature=0.8,
                                   title=title,clusterAlg="hc", distance="pearson",seed=1262118388.71279,plot="png")
  }
  # results = ConsensusClusterPlus(d,maxK=6,reps=50,pItem=0.8,pFeature=1,
  #                                + title=title,clusterAlg="hc",distance="pearson",seed=1262118388.71279,plot="png")
  dir.create(here::here("results","cluster_compares", title))
  
  #get item consensus for the clusters
  icl = calcICL(results,title=title,plot="png")
  icl[["clusterConsensus"]]
  write.table(icl[["clusterConsensus"]],here::here("results","cluster_compares",title,"ClusterConsensus_Assignments.tsv"), sep = "\t", col.names = TRUE, row.names = FALSE)
  
  #For every cluster tested, write the sample cluster assignments
  for(i in 2:7){#TODO is this a list?? change to lapply
    cluster_temp <- data.frame("Samples" = names(results[[i]][["consensusClass"]]), "ClusterID" = results[[i]][["consensusClass"]])
    cluster_temp$ClusterID <- paste0("Cluster", cluster_temp$ClusterID)
    write.table(cluster_temp,here::here("results","cluster_compares",title,paste0("Sample_Cluster_Assignments_",i,".tsv")), sep = "\t", col.names = TRUE, row.names = FALSE)
    
  }
  message("FINISHED", toupper(cluster_type),"CLUSTERING","\n")
  message("STARTING SILOUETTE CALCULATIONS FOR", toupper(cluster_type),"CLUSTERING","\n")
  silouette_scores(swept_away, title, "normal")
  message("COMPLETED SILOUETTE CALCULATIONS FOR", toupper(cluster_type),"CLUSTERING","\n")
}


  expr_dat <- readRDS(here::here("data","disease_expression.rds"))
  message("READ IN DATA\n")
  
  # 2.0 normal clustering ----
  library(ConsensusClusterPlus)
  
  #row-wise median centering
  swept_away <- sweep(expr_dat,1, apply(expr_dat,1,median,na.rm=T))
  message("NORMALIZED DATA\n")
  # Cluster with PAM, K Means and HC for comparisons
  
  Cluster("pam",swept_away, "PAM_Val")
  Cluster("hc",swept_away, "HC_Val")
  Cluster("km",swept_away, "KM_Val")
  
  
  
  
