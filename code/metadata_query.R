#####################################INFO#######################################
########### TITLE: FOR QUERYING GEO DATASETS 
########### PROJECT: SEPSIS META-ANALYSIS 
########### DESCRIPTION: SCRIPT TO QUERY DATASETS FROM GEO
########### AUTHOR: LESLIE SMITH 
########### DATE: JANUARY 15TH, 2022

library("optparse") 
library(GEOquery)
require(readr)  # for read_csv()
require(dplyr)  # for mutate()
require(tidyr)  # for unnest()
require(purrr)
library(here)

#input should be name of a series
option_list <- list(
  make_option(c("-d", "--dataset"), type="character", default=NULL, 
              help="Name of dataset to quantify", metavar="character")
); 

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
#Query dataset from GEO 
gse <- getGEO(study,GSEMatrix=TRUE)
gse_list <- gse[[1]]
#get phenotype data
metadata <- pData(phenoData(gse_list))
#get platform ID for mapping if needed
platform <- metadata$platform_id[1]
#write metadata to file to be combined in master metadata file
meta_row <- tibble::rownames_to_column(metadata,"Sample_ID")
write.table(meta_row,here::here("datasets","Metadata_Curation",paste0(study,"_Metadata.tsv")),sep = "\t",col.names = TRUE,row.names = FALSE)

#get expression data
expression <- exprs(gse_list)

#platform mapping example:
platform <- readr::read_tsv(here::here("datasets","GSE134347","GPL17586-45144.txt"), skip = 15)
# split gene column by " // " to get the symbol
plat_map <- platform |> mutate(temp = str_extract(gene_assignment, "(?<= // )[^/]+")) |> dplyr::select(ID, temp)
head(plat_map)
#get rid of the spaces
plat_cut <- plat_map |> mutate(SYMBOL = str_replace_all(temp, ' ',"")) |> dplyr::select(ID, SYMBOL)
head(plat_cut)

expression_dat <- tibble::rownames_to_column(as.data.frame(expression), "probeID")
expression_dat[1:5,1:5]

mapped <- merge(plat_cut, expression_dat, by.x = "ID", by.y = "probeID")
mapped[1:5,1:5]
mapped2 <- mapped[,-1]
mapped2[1:5,1:5]
mapped3 <- mapped2[mapped2$SYMBOL != "---",]
mapped3[1:5,1:5]
write.table(mapped3, here::here("datasets","GSE134347","GSE134347_expr_mapped.tsv"), sep = "\t", col.names = TRUE, row.names = FALSE)
