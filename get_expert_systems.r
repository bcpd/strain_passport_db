library(jsonlite)
library(tidyverse)
library(data.table)
library(readr)

# read MAG*.json files
mags <- list.files(pattern = "json", ignore.case = T, full.names = T)

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
write_csv(vfdb1, "vfdb-annotations.csv")

# 2.1. Create VFG-VFID-VFC mapping table from fasta file
library(tidysq)

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


# amrfinder

amrs <- dt %>%
  filter(!is.na(expert.amrfinder.gene))
         
#amr_hmms <- dt %>%
#  filter(!is.na(expert.amrfinder.gene) &
##           expert.amrfinder.method == "HMM")

#amr_blast <- dt %>%
#  filter(!is.na(expert.amrfinder.gene) &
 #          expert.amrfinder.method != "HMM")

amrcats_url <- "https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest/ReferenceGeneCatalog.txt"

amrcats <- fread(amrcats_url)

amrs <- amrs %>%
  rename(gene_family = expert.amrfinder.gene)

amrs_cats <- amrs %>%
  inner_join(amrcats, by = "gene_family") %>%
  distinct(locus, .keep_all = T)

#amrs_cats %>%
#  select(gene_family, class.y, subclass) %>% View()



#destfile <- getwd()

#download.file(amrcats_url, destfile)





