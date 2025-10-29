## reannotating with a new library:

library(EnsDb.Hsapiens.v86)
edb <- EnsDb.Hsapiens.v86

## GEO Datasets: Get to Symbol genes -> take the mean to combine duplicates -> QN in conda 
#GSE100159
GSE100159 <- readr::read_tsv(here::here("datasets","GSE100159","GSE100159_expression.txt"), col_names = T)
GSE100159[1:5,1:5]
#problems()
GSE100159 <- GSE100159[-48804,] #last row na
GSE100159_plat <- readr::read_tsv(here::here("datasets","GSE100159","GPL6884.txt"), col_names = T)
GSE100159_plat[1:5,1:5]
#Mc("row","Entrez", "ID", "Symbol", "enseml_gene_id"), skip = 1
GSE100159_plat <- GSE100159_plat %>% dplyr::select(ID, Symbol)
colnames(GSE100159_plat) <- c("ID", "SYMBOL")
GSE100159_plat[1:5,1:2]
GSE100159_map <- merge(x = GSE100159_plat, y = GSE100159, by.x = "ID", by.y = "ID_REF")
GSE100159_map[1:5,1:5]
GSE100159_map <- GSE100159_map[,-1] # get rid of probe ID
# make sure expression values are numeric
GSE100159_map[,2:length(GSE100159_map)] <- lapply(GSE100159_map[,2:length(GSE100159_map)], as.numeric)
GSE100159_map[1:5,1:5]
GSE100159.dat <- GSE100159_map %>% group_by(SYMBOL) %>% summarise_all(mean) # 25441 x 48 ??? 
GSE100159.dat[1:5,1:5]
dim(GSE100159.dat)
write.table(GSE100159.dat, file = here::here("datasets","GSE100159","Symbol_Expression.tsv"),sep = "\t", col.names = T, row.names = F)

#------------------------------- GSE106878 -------------------------------
#GSE106878 
GSE106878 <- readr::read_tsv(here::here("datasets","GSE106878","GSE106878_expression.txt"), col_names = T)
GSE106878 <- GSE106878[-29648,] #last row was na 
GSE106878_plat <- readr::read_tsv(here::here("datasets","GSE106878","GPL10295_trimmed.txt"))
GSE106878_plat[1:5,1:5]
GSE106878_plat <- GSE106878_plat %>% dplyr::select(ID,Symbol)
colnames(GSE106878_plat) <- c("ID","SYMBOL")
GSE106878_plat[1:5,1:2]
GSE106878_map <- merge(x = GSE106878_plat, y = GSE106878, by.x = "ID", by.y = "ID_REF")
GSE106878_map[1:5,1:5]
GSE106878_map <- GSE106878_map[,-1] #get rid of probe id
GSE106878.dat <- GSE106878_map %>% group_by(SYMBOL) %>% summarise_all(mean) #dim 18976    95
dim(GSE106878.dat)
write.table(GSE106878.dat, file = here::here("datasets","GSE106878","Symbol_Expression.tsv"),sep = "\t", col.names = T, row.names = F)

#------------------------------- GSE131761 -------------------------------
GSE131761 <- readr::read_tsv(here::here("datasets","GSE131761","GSE131761_trimmed.txt"), col_names = T)
GSE131761[1:5,1:5]
GSE131761 <- GSE131761[-34128,] #last row na
GSE131761_plat <- readr::read_tsv(here::here("datasets","GSE131761","GPL13497-9755_cut.txt"), col_names = T)
GSE131761_plat[1:5,1:5]
GSE131761_plat <- GSE131761_plat %>% dplyr::select(ID, GENE_SYMBOL)
colnames(GSE131761_plat) <- c("ID","SYMBOL")
GSE131761_map <- merge(x = GSE131761_plat, y = GSE131761, by.x = "ID", by.y = "ID_REF")
GSE131761_map[1:5,1:5]
GSE131761_map  <- GSE131761_map[,-1]#get rid of probe id
class(GSE131761_map$GSM3816544)
GSE131761.dat <- GSE131761_map %>% group_by(SYMBOL) %>% summarise_all(mean) #dim 21755   130
dim(GSE131761.dat)
write.table(GSE131761.dat, file = here::here("datasets","GSE131761","Symbol_Expression.tsv"),sep = "\t", col.names = T, row.names = F)

#------------------------------- GSE137340 -------------------------------
#### - platform mapping addition on May 8th 2024 - ####
dat <- readr::read_tsv(here::here("datasets","GSE137340","GSE137342-GPL10558_series_matrix.txt"))
platform15 <- readr::read_tsv(here::here("datasets","GSE137340","GPL10558-50081.txt"), skip = 30)
head(platform15)
colnames(platform15)
plat_map <- platform15 |> dplyr::select(ID,ILMN_Gene)
head(plat_map)
dat_mapped <- merge(plat_map,dat, by.x = "ID",by.y = "ID_REF")
dat_mapped[1:5,1:5]
dat_mapped <- dat_mapped[,-1] #drop ensembl
dat_mapped[1:5,1:5]
names(dat_mapped)[names(dat_mapped) == "ILMN_Gene"] <- "SYMBOL"
dat_mapped <- dat_mapped |> dplyr:: group_by(SYMBOL) |> dplyr::summarise_all(mean) 
write.table(dat_mapped, here::here("datasets","GSE137340","Symbol_Expression.tsv"),sep = "\t",col.names = TRUE, row.names = FALSE)

#------------------------------- GSE189400 -------------------------------
#note that this is the gene count subseries from GSE189403
#just mapping to SYMBOL here
dat <- readr::read_csv(here::here("datasets","GSE189400","GSE189400_earli_paxgene_counts_50k.csv"))
dat[1:5,1:5]
dat_mutated <- dplyr::mutate(dat, cut_ensembl = gsub("\\..*", "", dat$ENSEMBL)) #there are .123 numbers on the end of the genes
to_symbol <- AnnotationDbi::select(edb, keys = dat_mutated$cut_ensembl,
                                            keytype = "GENEID",
                                            columns = c("SYMBOL","GENEID"))
to_symbol[1:5,1:2]
master_dat <- merge(to_symbol, dat_mutated, by.x = "GENEID", by.y = "cut_ensembl")
master_dat[1:5,1:5]
dim(master_dat)
master_cut <- master_dat |> dplyr::select(-GENEID, -ENSEMBL)
master_cut[1:5,1:5]
to_write <- master_cut |> dplyr::group_by(SYMBOL) |> dplyr::summarise_all(mean)
#we had weird duplicate columns for some samples. I decided just to remove them, we lose 31 samples including non-sepsis samples.
df <- to_write[, !grepl("\\.\\.\\.", names(to_write))] 
dim(df)
write.table(df, here::here("datasets","GSE189400","Symbol_Expression.tsv"),sep = "\t", col.names = TRUE, row.names = FALSE)

#------------------------------- GSE236713 -------------------------------
platform <- readr::read_tsv(here::here("datasets","GSE236713","GPL17077-17467.txt"),skip = 16 ) 
head(platform)
dat <- readr::read_tsv(here::here("datasets","GSE236713","GSE236713_series_matrix_copy.txt"))
dat[1:5,1:5]
dat <- dat[dat$ID_REF != "!series_matrix_table_end",]
dim(dat)
map <- platform |> dplyr::select(ID,GENE_SYMBOL)
head(map)
dat_mapped <- merge(map, dat, by.x = "ID", by.y = "ID_REF")
dat_mapped[30700,1:5]
dat_cut <- dat_mapped[,-1]
dim(dat_cut)
dat_cut[1:5,1:5]
dat_cut <- as.data.frame(dat_cut)
dat_cut <-dat_cut[dat_cut$GENE_SYMBOL != "<NA>",]
names(dat_cut)[names(dat_cut) == "GENE_SYMBOL"] <- "SYMBOL"
dat_cut <- dat_cut[!(is.na(dat_cut$GSM7573957)),] #not sure why sometimes it keeps the unwanted rows and just turns the values to NA??
#it seems there are actually 4535 probe IDs that map to NA
dat <- dat_cut |> dplyr::group_by(SYMBOL) |> dplyr::summarise_all(mean) 

write.table(dat, here::here("datasets","GSE236713","Symbol_Expression.tsv"),sep = "\t", col.names = TRUE, row.names = FALSE)

#------------------------------- GSE236892 -------------------------------
#was already in Ensembl IDs
#already in ensembl genes
dat <- readr::read_csv(here::here("datasets","GSE236892","GSE236892_cnt_data.csv"))
names(dat)[names(dat) == "...1"] <- "ENSEMBL"
to_symbol <- AnnotationDbi::select(edb, keys = dat$ENSEMBL,
                    keytype = "GENEID",
                    columns = c("SYMBOL","GENEID"))
head(to_symbol)
master_dat <- merge(to_symbol, dat, by.x = "GENEID", by.y = "ENSEMBL")
master_dat[1:5,1:5]
dim(master_dat)
master_cut <- master_dat |> dplyr::select(-GENEID)
master_cut[1:5,1:5]
to_write <- master_cut |> dplyr::group_by(SYMBOL) |> dplyr::summarise_all(mean) 
write.table(to_write, here::here("datasets","GSE236892","Symbol_Expression.tsv"),sep = "\t",col.names = TRUE, row.names = FALSE)
#------------------------------- GSE32707 -------------------------------
GSE32707 <- readr::read_tsv(here::here("datasets","GSE32707","GSE32707_trimmed.txt"), col_names = T)
GSE32707 <- GSE32707[GSE32707$ID_REF != "!series_matrix_table_end",]
GSE32707_plat <- readr::read_tsv(here::here("datasets","GSE32707","GPL10558_trimmed.txt"), col_names = T)
GSE32707_plat[1:5,1:5]
GSE32707_t <- GSE32707_plat %>% dplyr::select(ID, ILMN_Gene)
GSE32707_t[1:5,1:2]
colnames(GSE32707_t) <- c("ID","SYMBOL")
GSE32707_t[1:5,1:2]
dim(GSE32707_t)
GSE32707_merg <- merge(x = GSE32707_t, y = GSE32707, by.x = "ID", by.y = "ID_REF")
GSE32707_merg[1:5,1:5]
GSE32707_merg <- GSE32707_merg[,-1] #get rid of the Probe ID
GSE32707.t <- GSE32707_merg %>% group_by(SYMBOL) %>% summarise_all(mean)
dim(GSE32707.t)
write.table(GSE32707.t, file = here::here("datasets","GSE32707","Symbol_Expression.tsv"),sep = "\t", col.names = T, row.names = F)

#------------------------------- GSE69063 -------------------------------
GSE69063 <- readr::read_tsv(here::here("datasets","GSE69063","GSE69063_trimmed.txt"), col_names = T)
GSE69063 <- GSE69063[-nrow(GSE69063),]
GSE69063_plat <- readr::read_tsv(here::here("datasets","GSE69063","GPL19983.txt"), col_names = T)
GPL19983_merg <- merge(x = GSE69063_plat, y = GSE69063, by.x = "ID", by.y = "ID_REF")
GPL19983_merg <- GPL19983_merg[,-1]
entrez_to_symbol <- AnnotationDbi::select(edb, keys = as.character(dat$ENTREZ_GENE_ID), keytype = "ENTREZID",columns = c("SYMBOL","ENTREZID"))
entrez_to_symbol[1:5,1:2]

GSE69063_map <- merge(x= entrez_to_symbol, y= GPL19983_merg, by.x = "ENTREZID", by.y = "ENTREZ_GENE_ID")
GSE69063_map <- GSE69063_map[,-1]
GSE69063.t <- GSE69063_map |> dplyr::group_by(SYMBOL) |> dplyr::summarise_all(mean)
GSE69063.t[1:5,1:5]
dim(GSE69063.t)
write.table(GSE69063.t, file = here::here("datasets","GSE69063","Symbol_Expression.tsv"),sep = "\t", col.names = T, row.names = F)

#------------------------------- GSE134347 ------------------------------- 
#using previosuly mapped expression matrix.

dat <- readr::read_tsv(here::here("datasets","GSE134347","GSE134347_expr_mapped.tsv"), col_names = T)

dat_mapped <- dat |> dplyr::group_by(SYMBOL) |> dplyr::summarise_all(mean)
write.table(dat_mapped, file = here::here("datasets","GSE134347","Symbol_Expression.tsv"),sep = "\t", col.names = T, row.names = F)


platform <- readr::read_tsv(here::here("datasets","GSE134347","GPL17586-45144.txt"), skip = 15)
platform_mut <- dplyr::mutate(platform, SYMBOL = get_gene(gene_assignment))

# ------------------------------- EMTAB1548 -------------------------------
# messy because of first line read in caused everything to be read in as a character 
EMTAB1548 <- readr::read_tsv(here::here("datasets","EMTAB1548","Exp_Sepsis_2014.txt"))
EMTAB1548[1:5,1:5]
EMTAB1548 <- EMTAB1548[-1,] #remove first row which contained reporter information

EMTAB1548_plat <- readr::read_tsv(here::here("datasets","EMTAB1548","A-MEXP-2183_trimmed.txt"),col_names = FALSE)
EMTAB1548_plat[1:5,1:11]
#map transcripts to ensembl genes
EMTAB1548_plat_trim <- EMTAB1548_plat |> dplyr::select(X5, X11)
head(EMTAB1548_plat_trim)
#If multiple transcripts - take the first one 
map2 <- EMTAB1548_plat_trim |> dplyr::mutate(trimmed = sapply(stringr::str_split(X11, ";"),`[`,1))
head(map2)

#get rid of NAs
map_final <- map2[!(is.na(map2$trimmed)),]
# get genes mappings
transcript_to_gene <- AnnotationDbi::select(edb, 
                             keys = map_final$trimmed,
                             keytype = "TXID",
                             columns = c("SYMBOL","GENEID", "TXID"))
head(transcript_to_gene)
# map probes to ensembl transcripts and genes
all_map <- merge(map_final, transcript_to_gene, by.x = "trimmed", by.y = "TXID") 
head(all_map)
all_map <- all_map[,-3] #we dont need the old transcript column
# merge with expression matrix
expr_dat <- merge(all_map, EMTAB1548, by.x = "X5", by.y = "Hybridization REF")
expr_dat[1:5,1:5]
# get rid of symbol column for now
expr <- expr_dat |> dplyr::select(-X5, -trimmed, -GENEID) # 27629 x 154
expr[1:5,1:5]
#group by gene and take average expression (this is same as what Refine.Bio does)
expr[,2:length(expr)] <- lapply(expr[,2:length(expr)], as.numeric)
expr_mod <- expr |> dplyr::group_by(SYMBOL) |> dplyr::summarise_all(mean) #dim 16608 x 154
dim(expr_mod)
write.table(expr_mod, here::here("datasets","EMTAB1548","Symbol_Expression.tsv"), sep = "\t", col.names = TRUE, row.names = FALSE)



#------------------------------- EMTAB4421 -------------------------------
EMTAB4421 <- readr::read_tsv(here::here("datasets","EMTAB4421","Davenport_sepsis_Jan2016_normalised_265.txt"), col_names = T)
EMTAB4421_plat <- readr::read_tsv(here::here("datasets","EMTAB4421","A-MEXP-2210_trimmed.txt"), col_names = F)
head(EMTAB4421_plat)
EMTAB4421_plat_trim <- EMTAB4421_plat %>% dplyr::select(X1, X5)
first <- merge(x = EMTAB4421_plat_trim, y = EMTAB4421, by.x = "X1", by.y = "Probe_ID")
first[1:5,1:5]
first <- first[,-1]#get rid of probes

EMTAB4421_metadata <- readr::read_tsv(here::here("datasets","EMTAB4421","metadata_EMTAB4421.tsv"), col_names = T)
EMTAB4421_metadata[1:5,1:5]
# get sample ID thats in the metadata and map to the ones in the expression matrix
mets <- EMTAB4421_metadata %>% dplyr::select(`Source Name`, `Comment[beadchiparray_id]`)
mets[1:5,1:2]
names(first) <- plyr::mapvalues(colnames(first), from = mets$`Comment[beadchiparray_id]`, to = mets$`Source Name`)
first[1:5,1:5]
colnames(first)[colnames(first) == "X5"] <- "SYMBOL"
#get the expression values to be numeric
first[,2:length(first)] <- lapply(first[,2:length(first)], as.numeric)
final <- first |> group_by(SYMBOL) |> summarise_all(mean)#dim 18428 x 266
write.table(final,here::here("datasets","EMTAB4421","Symbol_Expression.tsv"), sep = "\t", col.names = TRUE, row.names = FALSE)

# ------------------------------- GSE13015 -------------------------------
plat1 <- readr::read_tsv(here::here("datasets","GSE13015","GPL6106-11578.txt"))
plat1[1:5,1:5]
plat1_cut <- plat1 |> dplyr::select(ID, Symbol)
head(plat1_cut)
plat2 <- readr::read_tsv(here::here("datasets","GSE13015","GPL6947-13512.txt"))
colnames(plat2)
head(plat2$Symbol)
plat2_cut <- plat2 |> dplyr::select(ID, Symbol)
head(plat2_cut)
dat1 <- readr::read_tsv(here::here("datasets","GSE13015","GSE13015-GPL6106_trimmed.txt"))
dat1 <- dat1[dat1$ID_REF != "!series_matrix_table_end",]
dat2 <- readr::read_tsv(here::here("datasets","GSE13015","GSE13015-GPL6947_trimmed.txt"))
dat2 <- dat2[dat2$ID_REF != "!series_matrix_table_end",]

merge1 <- merge(plat1_cut, dat1, by.x = "ID", by.y = "ID_REF")
dim(merge1)
merge1[1:5,1:5]
merge1 <- merge1[,-1]
names(merge1)[names(merge1) == "Symbol"] <- "SYMBOL"
final1 <- merge1 |> group_by(SYMBOL) |> summarise_all(mean)
final1[1:5,1:5]
dim(final1)

merge2 <- merge(plat2_cut, dat2, by.x = "ID", by.y = "ID_REF")
dim(merge2)
merge2[1:5,1:5]
merge2 <- merge2[,-1]
names(merge2)[names(merge2) == "Symbol"] <- "SYMBOL"
final2 <- merge2 |> group_by(SYMBOL) |> summarise_all(mean)
final2[1:5,1:5]
dim(final2)

all <- merge(final1, final2, by = "SYMBOL")
dim(all)
write.table(all,here::here("datasets","GSE13015","Symbol_Expression.tsv"), sep = "\t", col.names = TRUE, row.names = FALSE)


