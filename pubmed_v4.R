#rm(list=ls())


# Config
# Change this to point to the location of the gtdbtk taxonomic classification

project_code_file <- "Project_code.txt"
MAGS_taxonomy_file <- "gtdbtk.bac120.summary.tsv"
api_key <- ""

# Required libraries
library(tidyverse)
library(readr)
#library(plyr)
library("easyPubMed")


# Reading the project ID
# This is important to make each MAGID unique as the output for multiple projects can be the same e.g MAG??

project_id <- read.delim(file = project_code_file, header = FALSE)
project_id <- as.character(project_id$V1)



# Read MAG taxonomy file to create list of organisms names to search

mags <- readr::read_tsv(MAGS_taxonomy_file)
mag_sp <- mags$classification %>%
  str_split(";", simplify = T) %>%
  as_tibble() %>%
  pull("V7") %>%
  str_replace_all(c("s__|_A|_C|_E|Imtechella halotolerans"), "") %>%
  str_replace("Glutamicibacter sp004320535", "Arthrobacter sp. S41") %>%
  na_if("") %>%
  as_tibble() %>%
  #filter(!is.na(value)) %>%
  pull("value")

mag_sp[7]  = "Arthrobacter S41"
mag_sp2 = data.frame(MAGID = paste0(project_id, "_", mags$user_genome),
                     Species = mag_sp)


# Debuging block
#query_year=2020
#the_species = "Rhodobacter sphaeroides"
#Genus = "Rhodobacter"
#Sp = "sphaeroides"
#counter = 5

counter = 1
while(counter <= nrow(mag_sp2)) {
    the_species = as.character(mag_sp2$Species[counter])
    if (is.na(the_species)){
      # Do nothing
      } else {
        
        print(paste("Working on", the_species))
        
        Genus = data.frame(str_split(the_species,  pattern = " "))[1,]
        Sp = data.frame(str_split(the_species,  pattern = " "))[2,]
        
        for (query_year in seq(from = 2001, to = 2022, by =1 )){
          Sys.sleep(4)
          print(paste("Year:", query_year))
          
          new_query = paste0( Genus, "[TIAB] AND ",
                              Sp, "[TIAB] AND ",
                              query_year, "[DP]")
          new_query <- get_pubmed_ids(new_query,
                                      api_key = api_key
                                          )
          
          # Check number of results and download data
          print(paste(new_query$Count, "records found"))
          
          if (new_query$Count ==0) {
            # Do nothing
          } else {
            fetched_data <- fetch_pubmed_data(new_query,
                                              encoding = "ASCII", 
                                              retmax = 4999)
            new_PM_df <- table_articles_byAuth(pubmed_data = fetched_data, 
                                               included_authors = "first", 
                                               max_chars = 1000, 
                                               encoding = "ASCII")
            new_PM_df$abstract = gsub('[\n\r]',
                                      ' ',
                                      x = new_PM_df$abstract)
            new_PM_df$address = NULL 
            new_PM_df$email = NULL
            
            fileout_fp = as.character(paste0(mag_sp2$MAGID[counter],
                                            "_",
                                            query_year,
                                            "_pubmed.txt"))
            write.table(x=new_PM_df,
                      sep = "\t",
                      quote = FALSE,
                        row.names = FALSE, 
                          file=fileout_fp)
          }
        }
      }
    counter = counter + 1
}



## Merge_files
## Get list of MAGS that do not have na in their species

MAGID_list = subset(mag_sp2, !is.na(mag_sp2$Species))$MAGID

# Debug block
#MAGID = as.character(MAGID_list[6])
#MAGID

suppressWarnings(rm(dataset))
for (MAGID in MAGID_list) {
  MAGID = as.character(MAGID)
  print(paste("Working on ", MAGID))

  file_list <- list.files(path = getwd(),
                          pattern = paste0(MAGID,"_")
                          )
  # Debug
  #file_list
  suppressWarnings(rm(dataset))
  for (file in file_list){
    print(file)

        # if the merged dataset doesn't exist, create it
    if (!exists("dataset")){
      dataset <- read.table(file, header=TRUE, sep="\t", quote = "",fill=TRUE)
    }
    
    # if the merged dataset does exist, append to it
    if (exists("dataset")){
      temp_dataset <-read.table(file, header=TRUE, sep="\t", quote = "", 
                                fill=TRUE)
      dataset<-rbind(dataset, temp_dataset)
      rm(temp_dataset)
    }
    
    clean_dataset <- dataset %>%
      distinct(pmid, .keep_all = TRUE)
    clean_dataset$address = NULL
    clean_dataset$email = NULL
    fileout_fp2 = as.character(paste0(MAGID,
                                     "_all_pubmed.txt"))
    clean_dataset$abstract = gsub('[\n\r]',' ',x = clean_dataset$abstract)
    
    write.table(x=clean_dataset,
                sep = "\t",
                quote = FALSE,
                row.names = FALSE, 
                file=fileout_fp2)
    
  }
}


# Yearly files are removed manually




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
# [1] easyPubMed_2.13     medrxivr_0.0.5.9000 forcats_0.5.0       stringr_1.4.1       dplyr_1.0.2        
# [6] purrr_0.3.4         readr_1.3.1         tidyr_1.1.2         tibble_3.0.3        ggplot2_3.3.5      
# [11] tidyverse_1.3.0    
# 
# loaded via a namespace (and not attached):
# [1] Rcpp_1.0.9        pillar_1.4.6      compiler_3.5.3    cellranger_1.1.0  dbplyr_1.4.4      prettyunits_1.1.1
# [7] progress_1.2.2    tools_3.5.3       jsonlite_1.8.0    lubridate_1.8.0   lifecycle_1.0.1   gtable_0.3.0     
# [13] pkgconfig_2.0.3   rlang_0.4.11      reprex_0.3.0      cli_3.2.0         rstudioapi_0.13   DBI_1.1.0        
# [19] curl_4.3.2        haven_2.3.1       withr_2.5.0       xml2_1.3.2        httr_1.4.4        fs_1.5.0         
# [25] generics_0.1.3    vctrs_0.3.4       hms_0.5.3         grid_3.5.3        tidyselect_1.1.0  glue_1.6.2       
# [31] R6_2.5.1          fansi_1.0.3       readxl_1.3.1      modelr_0.1.8      blob_1.2.1        magrittr_2.0.3   
# [37] backports_1.1.10  scales_1.1.1      ellipsis_0.3.2    rvest_0.3.6       assertthat_0.2.1  colorspace_2.0-3 
# [43] utf8_1.2.2        stringi_1.5.3     munsell_0.5.0     broom_0.8.0       crayon_1.5.1     
# > 