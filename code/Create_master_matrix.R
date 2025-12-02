# -- Read in all datasets and metadata and put them into 1 common expression matrix and 1 common metadata file. -- 
# -- This script gets us a master matrix with ALL samples from included datasets - there is no disease exclusion at this point -- 


# 1.0 Load All Expression Matrices ---- 

## define datasets to read in
studies <- list.files(here::here("datasets"))
## we dont include these so remove them (SRP have been replaced with our own seqeuencing):
#EMTAB1548, EMTAB4421, GSE10474, GSE106878 have been eliminated due to gene inclusion
studies <- studies[!(studies %in% c("README.txt", "salmon_quant","Metadata_Curation",
                                    "SRP049820","TEMP.txt", "additional_genes.txt", "Metadata.xlsx", 
                                    "EMTAB1548","EMTAB4421", "GSE10474", "GSE106878","SRP132709","Metadata.csv", "GSE186054", "GSE28750",
                                    "Spearman_Sample_Cor.tsv"
                                    ))]
## we only want study files here
studies <- studies[!grepl("Master*", studies)]
studies <- studies[!grepl("Disease*", studies)]
studies
length(studies) ##should now be 29
studies_paths <- here::here("datasets", studies,"QN_Expression.tsv")
##reading in all expression matrices into list
study_list <- lapply(studies_paths, readr::read_tsv) 


# 2.0 Create Batch IDs ----

#create batch IDs to use later on in Combat
batchIDs <- unlist(sapply(seq_along(studies), function(i) {rep(studies[i],ncol(study_list[[i]])-1)} ))
head(batchIDs)
length(batchIDs) ## should match the number of samples in the expresstion matrix
master <- study_list |> purrr::reduce(dplyr::inner_join) #10150  3713
master_row <- tibble::column_to_rownames(master, "SYMBOL") 
dim(master_row)
names(batchIDs) <- colnames(master_row)
head(batchIDs)
batchIDs <- batchIDs[names(batchIDs) != "SR191"] # get rid of this sample becuase of low quality 
master_row <- master_row |>  dplyr::select(-SR191) 

## * -- save RDS -- *
saveRDS(batchIDs, here::here("data","batchIDs_allSamples.rds"))
saveRDS(master_row, here::here("data","raw_all_samples_expression.rds"))
