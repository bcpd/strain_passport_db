
# Config
# Change this to point to the location of the gtdbtk taxonomic classification

project_code_file <- "Project_code.txt"
bakta_folder = "MAGS_bakta"

#MAGS_taxonomy_file <- "gtdbtk.bac120.summary.tsv"
#api_key <- "0cb41d0e0c2136aff7e0478bc4182be1cb09"



# Reading the project ID
# This is important to make each MAGID unique as the output for multiple projects can be the same e.g MAG??

project_id <- read.delim(file = project_code_file, header = FALSE)
project_id <- as.character(project_id$V1)


library(jsonlite)
library(tidyverse)
library(data.table)
library(readr)
library(tidysq)


# read MAG*.json files
mags <- list.files(path = bakta_folder,
                   pattern = "json",
                   ignore.case = T,
                   full.names = T)

maglist <- setNames(lapply(mags, read_json, simplifyVector = T, flatten = T),
         tools::file_path_sans_ext(basename(mags)))

ef <- enframe(maglist)

efwid <- ef %>%
  unnest_auto(value) %>%
  select(name, features)

setattr(efwid$features, 'names', efwid$name)

dtfull <- data.table::rbindlist(efwid$features, fill = TRUE,
                                use.names = T, idcol = "SRA")

# transform list-columns into character vectors

dtfullist <- dtfull %>%
  select_if(is.list) %>%
  map_depth(2, paste0, collapse = ";") %>%
  map_dfr(unlist) %>%
  na_if("") %>%
  select(-pfams)

dtfull_notlist <- dtfull %>%
  select_if(negate(is.list))

dt <- tibble(dtfull_notlist, dtfullist)

# split table by expert system
# VFDB
# 1. create "expert.expert_proteins.source" VFDB table

vfdb <- dt %>%
  filter(expert.expert_proteins.source == "VFDB")

# 2. separate VFDB table into VFG and VFC elements
vfdb1 <- vfdb %>%
  separate(expert.expert_proteins.db_xrefs, into = c("VFG", "VFC"),
           remove = F, sep = ";") %>%
  mutate(
    VFG = str_replace(VFG, "VFDB:", ""),
    VFC = str_replace(VFC, "VFDB:", "")
      )
vfdb1$MAGID = paste0(project_id, "_", vfdb1$SRA)
write_csv(vfdb1, "vfdb-annotations.csv")



# Downloaded on 2022.10.06 from here both fasta and excel files from:
# http://www.mgc.ac.cn/VFs/download.htm

# 2.1. Create VFG-VFID-VFC mapping table from fasta file
#library(tidysq)
#update.packages("Rcpp")

vfdb_fas_headers <- read_fasta("VFDB_setA_nt.fas")

vfcats <- readxl::read_excel("VFs.xls", skip = 1, col_names = T)

vfcats <- vfcats %>%
  unite("VFSubcategory", VF_FullName:VF_Name,
        na.rm = TRUE, remove = FALSE)

vfdb_fas_headers <- vfdb_fas_headers %>%
  select(-sq) %>%
  mutate(
    VFG = str_extract(vfdb_fas_headers$name, "VFG\\d+"),
    VFID = str_extract(vfdb_fas_headers$name, "VF\\d{4}"),
    VFCID = str_extract(vfdb_fas_headers$name, "VFC\\d{4}"),
)

# join tables
vfdb_map <- vfdb_fas_headers %>%
  left_join(vfcats, by = c("VFCID", "VFID")) %>%
  select(-name, -Reference)

write_csv(vfdb_map, "vfdb_map.csv", na = "")



# Data from NCBI Antimicrobial resistance database
# Amrfinder

amrs <- dt %>%
  filter(!is.na(expert.amrfinder.gene))


amrcats_url <- "https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest/ReferenceGeneCatalog.txt"

amrcats <- fread(amrcats_url)

amrs <- amrs %>%
  rename(gene_family = expert.amrfinder.gene)

amrs_cats <- amrs %>%
  inner_join(amrcats, by = "gene_family") %>%
  distinct(locus, .keep_all = T)

amrs_cats$MAGID = paste0(project_id, "_", amrs_cats$SRA)
amrs_cats %>%
  select(MAGID, SRA, locus, pscc.product, gene_family, class.y, subclass) %>%
  write_csv("armfinder_annotation.csv", na="")

amrs_cats %>%
  select(gene_family, class.y, subclass) %>%
  write_csv("amrfinderplus_map.csv", na = "")


# sessionInfo()
# R version 4.0.3 (2020-10-10)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 19044)
# 
# Matrix products: default
# 
# Random number generation:
# RNG:     Mersenne-Twister 
# Normal:  Inversion 
# Sample:  Rounding 
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
# [1] data.table_1.13.6 forcats_0.5.1     stringr_1.4.0     dplyr_1.0.3       purrr_0.3.4       readr_1.4.0      
# [7] tidyr_1.1.2       tibble_3.0.5      ggplot2_3.3.3     tidyverse_1.3.0   jsonlite_1.7.2    tidysq_1.1.3     
# 
# loaded via a namespace (and not attached):
# [1] Rcpp_1.0.9        cellranger_1.1.0  pillar_1.4.7      compiler_4.0.3    dbplyr_2.1.0      tools_4.0.3      
# [7] lubridate_1.7.9.2 lifecycle_0.2.0   checkmate_2.0.0   gtable_0.3.0      pkgconfig_2.0.3   rlang_0.4.10     
# [13] reprex_1.0.0      cli_2.2.0         DBI_1.1.1         rstudioapi_0.13   curl_4.3          haven_2.3.1      
# [19] withr_2.4.1       xml2_1.3.2        httr_1.4.2        fs_1.5.0          hms_1.0.0         generics_0.1.0   
# [25] vctrs_0.3.6       grid_4.0.3        tidyselect_1.1.0  glue_1.4.2        R6_2.5.0          fansi_0.4.2      
# [31] readxl_1.3.1      modelr_0.1.8      magrittr_2.0.1    backports_1.2.1   scales_1.1.1      ellipsis_0.3.1   
# [37] rvest_0.3.6       assertthat_0.2.1  colorspace_2.0-0  stringi_1.5.3     munsell_0.5.0     broom_0.7.4      
# [43] crayon_1.3.4     
