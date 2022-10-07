# Config
# Change this to point to the location of the gtdbtk taxonomic classification

project_code_file <- "Project_code.txt"
taxonomy_file = "gtdbtk.bac120.summary.tsv"


# Reading the project ID
# This is important to make each MAGID unique as the output for multiple projects can be the same e.g MAG??

project_id <- read.delim(file = project_code_file, header = FALSE)
project_id <- as.character(project_id$V1)


# Required libraries

library(medrxivr)
library(tidyverse)


# load MAG taxonomies
mags <- readr::read_tsv(taxonomy_file)

# select species and remove I. halotolerans
mag_sp <- mags$classification %>%
  str_split(";", simplify = T) %>%
  as_tibble() %>%
  pull("V7") %>%
  str_replace_all(c("s__|_A|_C|_E|Imtechella halotolerans"), "") %>%
  str_replace("Glutamicibacter sp004320535", "Arthrobacter sp. S41") %>%
  na_if("") %>% as_tibble() %>% filter(!is.na(value)) %>% pull("value")



#change c. sphaer. to rhodo spha
all_sps <- str_replace(mag_sp, "Cereibacter sphaeroides",
                       "Rhodobacter sphaeroides")

# biorxiv database

#preprint_data_biorxiv <- mx_api_content(server = "biorxiv")
#saveRDS(preprint_data_biorxiv, "biorxiv_data_20221006")
preprint_data_biorxiv <- readRDS("biorxiv_data_20221006")


biorxiv_res <- all_sps %>%
  purrr::set_names() %>% #if not used, only stores index number
  purrr::map_df(function(x) mx_search(data = preprint_data_biorxiv,
                                     query = x, auto_caps = T),
                  .id = "species")

#replace "NA" with (true) NA
biorxiv_res <- biorxiv_res %>%
  mutate(across(where(is.character), ~na_if(., "NA")))

# unique(biorxiv_res$species) 15/22 species with results

## save biorxiv reults
write_csv(biorxiv_res, "biorxiv.csv")


#cs <- mx_search(data = preprint_data_biorxiv,
#          query = "Cereibacter sphaeroides", auto_caps = T)
#
#write_csv(cs, "biorxiv.csv")


# > sessionInfo()
# R version 4.1.3 (2022-03-10)
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
# [1] forcats_0.5.1   stringr_1.4.0   dplyr_1.0.8     purrr_0.3.4     readr_2.1.2     tidyr_1.2.0     tibble_3.1.6   
# [8] ggplot2_3.3.5   tidyverse_1.3.2 medrxivr_0.0.5 
# 
# loaded via a namespace (and not attached):
# [1] cellranger_1.1.0    pillar_1.8.1        compiler_4.1.3      dbplyr_2.1.1        prettyunits_1.1.1  
# [6] progress_1.2.2      tools_4.1.3         bit_4.0.4           lubridate_1.8.0     jsonlite_1.8.0     
# [11] googledrive_2.0.0   lifecycle_1.0.2     gargle_1.2.0        gtable_0.3.0        pkgconfig_2.0.3    
# [16] rlang_1.0.5         reprex_2.0.1        DBI_1.1.2           cli_3.2.0           rstudioapi_0.13    
# [21] curl_4.3.2          parallel_4.1.3      haven_2.4.3         xml2_1.3.3          withr_2.5.0        
# [26] httr_1.4.4          fs_1.5.2            generics_0.1.3      vctrs_0.4.0         hms_1.1.2          
# [31] bit64_4.0.5         googlesheets4_1.0.0 grid_4.1.3          tidyselect_1.1.2    glue_1.6.2         
# [36] data.table_1.14.2   R6_2.5.1            fansi_1.0.3         readxl_1.4.0        vroom_1.5.7        
# [41] modelr_0.1.8        tzdb_0.3.0          magrittr_2.0.3      backports_1.4.1     scales_1.1.1       
# [46] ellipsis_0.3.2      rvest_1.0.2         assertthat_0.2.1    colorspace_2.0-3    utf8_1.2.2         
# [51] stringi_1.7.6       munsell_0.5.0       broom_0.7.12        crayon_1.5.1     

