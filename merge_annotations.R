# Merge annotations
# This script processes annotation tables from bakta and dram

# You need to change this to process your file

project_code_file <- "Project_code.txt"
MAGS_bakta_folder <- "MAGS_bakta"
MAGS_DRAM_annotation_file <- "MAGS_DRAM/annotations/annotations.tsv"


## Required libraries
library(tidyverse)


## Reading the project ID
## This is important to make each MAGID unique as the output for
## multiple projects can be the same e.g MAG01 from project A and
## MAG01 from project B

project_id <- read.delim(file = project_code_file, header = FALSE)
project_id <- as.character(project_id$V1)


## This block reads the anotations from DRAM
## It also separates them into subtables just if needed

dram <- read.delim2(MAGS_DRAM_annotation_file)
names(dram)[1] <- "GeneID"
dram$genome_id <- paste0(project_id, "_", dram$fasta)
dram$fasta <- dram$genome_id
dram$genome_id <- NULL
names(dram)[2] <- "genome_id"
dram$heme_regulatory_motif_count <- NULL


dram_ko <- dram %>%
  select(genome_id, GeneID, ko_id) %>%
  filter(ko_id != "")

dram_kegg <- dram %>%
  select(genome_id, GeneID, kegg_hit) %>%
  filter(kegg_hit != "")

dram_pfam <- dram %>%
  select(genome_id, GeneID, pfam_hits) %>%
  filter(pfam_hits != "")

dram_cazy <- dram %>%
  select(genome_id, GeneID, cazy_id, cazy_hits) %>%
  filter(cazy_id != "") %>%
  filter(cazy_hits != "") %>%
  pivot_longer(cols = cazy_id:cazy_hits)


dram_peptidase <- dram %>%
  select(
    genome_id, GeneID, peptidase_id, peptidase_family, peptidase_hit,
    peptidase_RBH, peptidase_identity,
    peptidase_bitScore, peptidase_eVal
  ) %>%
  filter(peptidase_id != "") %>%
  pivot_longer(cols = peptidase_id:peptidase_eVal)


## Processing of bakta annotations
## Reads all the invidual MAG annotation tables, also separates the
## results from multiple databases into individual databases columns

baktas <- list.files(path = MAGS_bakta_folder,
                     pattern = "MAG[0-9]+.tsv", full.names = T)
bakta_all <- as.data.frame(matrix(ncol = 10))
colnames(bakta_all)[10] <- "MAGID"

for (bakta_file in baktas) {
  temp_bakta <- read.delim(file = bakta_file, header = F, comment.char = "#")
  bakta_file <- gsub(bakta_file, pattern = "[a-zA-Z0-9_/]+/", replacement = "")
  temp_bakta$MAGID <- gsub(x = bakta_file, pattern = ".tsv", replacement = "")
  bakta_all <- rbind(bakta_all, temp_bakta)
}


colnames(bakta_all) <- c(
  "Sequence_Id", "Type",
  "start_position", "end_position", "Strand",
  "Locus_Tag", "Gene", "Product", "DbXrefs",
  "MAGID"
)

bakta <- bakta_all[2:nrow(bakta_all), c(
  "MAGID", "Sequence_Id", "Type",
  "start_position", "end_position",
  "Strand", "Locus_Tag", "Gene",
  "Product", "DbXrefs"
)]

bakta$COG <- NA
bakta$KEGG <- NA
bakta$SO <- NA
bakta$UniParc <- NA
bakta$UniRef <- NA
bakta$RefSeq <- NA


counter <- 1
while (counter <= nrow(bakta)) {
  bakta_test <- bakta$DbXrefs[counter]
  bakta_test2 <- strsplit(x = bakta_test, split = ", ")
  bakta_test3 <- as.data.frame(bakta_test2)
  names(bakta_test3) <- "Variable"
  bakta_test4 <- separate(bakta_test3,
                          col = Variable, sep = ":", into = c("DB", "Value"))

  COG_set <- as.list(subset(bakta_test4, bakta_test4$DB == "COG")$Value)
  KEGG_set <- as.list(subset(bakta_test4, bakta_test4$DB == "KEGG")$Value)
  SO_set <- as.list(subset(bakta_test4, bakta_test4$DB == "SO")$Value)
  UniParc_set <- as.list(subset(bakta_test4, bakta_test4$DB == "UniParc")$Value)
  UniRef_set <- as.list(subset(bakta_test4, bakta_test4$DB == "UniRef")$Value)
  RefSeq_set <- as.list(subset(bakta_test4, bakta_test4$DB == "RefSeq")$Value)

  if (length(COG_set) == 0) {
    COG_set <- NA
  }
  if (length(KEGG_set) == 0) {
    KEGG_set <- NA
  }
  if (length(SO_set) == 0) {
    SO_set <- NA
  }
  if (length(UniParc_set) == 0) {
    UniParc_set <- NA
  }
  if (length(UniRef_set) == 0) {
    UniRef_set <- NA
  }
  if (length(RefSeq_set) == 0) {
    RefSeq_set <- NA
  }

  
  COG_set <- toString(COG_set)
  KEGG_set <- toString(KEGG_set)
  SO_set <- toString(SO_set)
  UniParc_set <- toString(UniParc_set)
  UniRef_set <- toString(UniRef_set)
  RefSeq_set <- toString(RefSeq_set)

  bakta$COG[counter] <- COG_set
  bakta$KEGG[counter] <- KEGG_set
  bakta$SO[counter] <- SO_set
  bakta$UniParc[counter] <- UniParc_set
  bakta$UniRef[counter] <- UniRef_set
  bakta$RefSeq[counter] <- RefSeq_set

  counter <- counter + 1
}



## Join bakta and DRAM tables

# Joins the bakta and DRAM annotation tables based on the GeneID which is a combination
# Of the MAGID and the Locus tag.

bakta$GeneID <- paste0(bakta$MAGID, "_", bakta$Locus_Tag)
all_annotations <- merge(bakta, dram, by = "GeneID")
write.table(
  x = all_annotations,
  file = "merged_annotations.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# Session infor
# If you have problems running the script you can try using the configuration
# Used to develop it


# sessionInfo()
# R version 3.5.3 (2019-03-11)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 19044)
# 
# Matrix products: default
# 
# locale:
# [1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252   
# [3] LC_MONETARY=English_United States.1252 LC_NUMERIC=C                          
# [5] LC_TIME=English_United States.1252    
# 
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] forcats_0.5.0   stringr_1.4.1   dplyr_1.0.2     purrr_0.3.4     readr_1.3.1     tidyr_1.1.2     tibble_3.0.3   
# [8] ggplot2_3.3.5   tidyverse_1.3.0
# 
# loaded via a namespace (and not attached):
#   [1] Rcpp_1.0.9       pillar_1.4.6     compiler_3.5.3   cellranger_1.1.0 dbplyr_1.4.4     tools_3.5.3     
# [7] jsonlite_1.8.0   lubridate_1.8.0  lifecycle_1.0.1  gtable_0.3.0     pkgconfig_2.0.3  rlang_0.4.11    
# [13] reprex_0.3.0     cli_3.2.0        rstudioapi_0.13  DBI_1.1.0        haven_2.3.1      withr_2.5.0     
# [19] xml2_1.3.2       httr_1.4.4       fs_1.5.0         generics_0.1.3   vctrs_0.3.4      hms_0.5.3       
# [25] grid_3.5.3       tidyselect_1.1.0 glue_1.6.2       R6_2.5.1         readxl_1.3.1     modelr_0.1.8    
# [31] blob_1.2.1       magrittr_2.0.3   backports_1.1.10 scales_1.1.1     ellipsis_0.3.2   rvest_0.3.6     
# [37] assertthat_0.2.1 colorspace_2.0-3 stringi_1.5.3    munsell_0.5.0    broom_0.8.0      crayon_1.5.1  
