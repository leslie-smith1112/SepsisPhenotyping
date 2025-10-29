## Sample mapping for matrices rerun on SRA but have GEO metadata
# in each datset here - I move the previously written QN_Expression.tsv file to QN_Expression_pre_sample_mapping.tsv
# This was done to maintain final matrices of each dataset as QN_Expression.tsv.

# These datasets will be mapped so that samples match the metadata. Probably not the best way to do this,
#Trying not to alter the metadata at all because when we start editing it we have 
# descrepencies between Hipergator version and issues version tracking between local copy. 

# index 13: SRR12291410: GSE154918 -- -- 
#   index 14: SRR16193818: GSE185263 -- -- 
#   index 15: SRR16465804: GSE186054 -- -- 
#   index 16: EARLI_10447 -- this one is correct in the title column -- GSE236892 -- 
# ^^ same here: GSE189400
# index 17: SRR17883537: GSE196117 --
#   index 18: SRR18547537: GSE199816 --
#   index 19: SRR21048668: GSE211210 --
#   index 20: SRR22106730: GSE216902 --
#   index 21: SRR23016072: GSE222393 --
#   index 22: SRR24524122: GSE232404 --
#   index 23: SRR24633016: GSE232753 --
#   index25: EARLI_10447: EARLI_10447 -- CHECK WITH KILEY AND FAHEEM
# index 30: SRR1652895: GSE63311


## NOTE: sadly not all SRA samples have same column name, 
datasets <- c("GSE154918","GSE185263","GSE186054","GSE196117", "GSE199816", "GSE211210", "GSE216902",
              "GSE222393", "GSE232404", "GSE232753", "GSE63311")

# for every dataset that needs to be mapped:
for(i in 1:length(datasets)){
  current <- datasets[i]
  message(current)
  #read in datasets:
  SRA_meta <- readr::read_csv(here::here("datasets",current,paste0(current,"_SRA_Metadata.csv")))
  expression <- readr::read_tsv(here::here("datasets",current,"QN_Expression_pre_sample_mapping.tsv"))
  
  #not all GEO sample IDs have the same column name in the SRA map, so
  # first we change the column name to be the same 
  colnames(SRA_meta) <- sapply(colnames(SRA_meta), function(x) {
    ifelse(substr(SRA_meta[1,x], 1, 3) == "GSM", "GEO_Accession", x)
  })
  
  colnames(expression) <- plyr::mapvalues(colnames(expression), from = SRA_meta$Run, to = SRA_meta$GEO_Accession)
  expression[1:5,1:5]
  write.table(expression, here::here("datasets",current,"QN_Expression.tsv"), sep = "\t", col.names = TRUE, row.names = FALSE)
}



current <- "GSE236892"
## these 2 are different: "GSE236892","GSE189400",
metadata <- readxl::read_xlsx(here::here("datasets","Metadata.xlsx"))
expression <- readr::read_tsv(here::here("datasets",current,"QN_Expression_pre_sample_mapping.tsv"))
the_map <- data.frame("EARLI_ID" = metadata$Title[metadata$Dataset == current], "GEO_Accession" = metadata$`GEO Accession`[metadata$Dataset == current])
#get rid of the suffixes:
the_map$EARLI_ID <- sub("_(Hyper|Hypo)$", "", the_map$EARLI_ID)
#map the GEPO expression
colnames(expression) <- plyr::mapvalues(colnames(expression), from = the_map$EARLI_ID, to = the_map$GEO_Accession)
expression[1:5,1:5]
write.table(expression, here::here("datasets",current,"QN_Expression.tsv"), sep = "\t", col.names = TRUE, row.names = FALSE)

## Similar for GSE189400
current <- "GSE189400"
expression <- readr::read_tsv(here::here("datasets",current,"QN_Expression_pre_sample_mapping.tsv"))
the_map <- data.frame("EARLI_ID" = metadata$Title[metadata$Dataset == current], "GEO_Accession" = metadata$`GEO Accession`[metadata$Dataset == current])
head(the_map)
#get rid of the suffixes:
the_map$EARLI_ID <- sub("_(SepsisBSI|Sepsisperiph|No-Sepsis)$", "", the_map$EARLI_ID)
table(the_map$EARLI_ID)
head(the_map)
#map the GEPO expression
colnames(expression) <- plyr::mapvalues(colnames(expression), from = the_map$EARLI_ID, to = the_map$GEO_Accession)
expression[1:5,1:5]
write.table(expression, here::here("datasets",current,"QN_Expression.tsv"), sep = "\t", col.names = TRUE, row.names = FALSE)





