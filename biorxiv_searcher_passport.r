install.packages("medrxivr")
library(medrxivr)
library(tidyverse)



# load MAG taxonomies
mags <- readr::read_tsv("Genomes_files/gtdbtk.bac120.summary.txt")

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

preprint_data_biorxiv <- mx_api_content(server = "biorxiv")
saveRDS(preprint_data_biorxiv, "biorxiv_data")

preprint_data_biorxiv <- readRDS("biorxiv_data")


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


cs <- mx_search(data = preprint_data_biorxiv,
          query = "Cereibacter sphaeroides", auto_caps = T)
  
 
