# GSE10474 -- QN_ensembl_expression
# GSE13015 -- good to go
# GSE28750 -- QN_ensembl_expression
# GSE66890
# GSE65682
# QN_ensembl_expression.tsv
# GSE74224
# GSE57065
# GSE95233
# GSE33118

library("optparse") 
library(tximport)
#library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ensembldb)
library(here)
library(tidyr)
##NOTE: You must change the database you used depending on the transcriptome version you used to annotate your samples in aligner
library(EnsDb.Hsapiens.v86)
edb <- EnsDb.Hsapiens.v86

option_list <- list(
  make_option(c("-d", "--dataset"), type="character", default=NULL, 
              help="Name of dataset to quantify", metavar="character")
); 

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
current_dat <- opt$dataset
#QN_ensembl_expression.tsv is the QN'd matrix from Refine.Bio
dat <- readr::read_tsv(here::here("datasets",current_dat,"QN_ensembl_expression.tsv"))
dat[1:5,1:5]

transcript_to_gene <- select(edb, 
                             keys = dat$Gene,
                             keytype = "GENEID",
                             columns = c("SYMBOL","GENEID"))
head(transcript_to_gene)
mapped <- merge(transcript_to_gene, dat, by.x = "GENEID", by.y = "Gene")
mapped[1:5,1:5]
mapped <- mapped[,-1]
expression <- mapped |> dplyr::group_by(SYMBOL) |> dplyr::summarise_all(mean)
dim(expression)
write.table(expression,here("datasets",current_dat,"QN_Expression.tsv"),sep = "\t", col.names = TRUE,row.names=FALSE)

