# Patent search


# We searched for patents using the patentsview package in R.
# This packages allow for search in the US patent system using its
# dedicated interface
# We looked for patents that included both the genus and species terms
# in the patent text.
# Example:
# The following query searched for patents with both
# "Cereibacter" and "sphaeroides" in the text:
# {"_and":[{"_text_all":{"patent_abstract":"Cereibacter"}},{"_text_all":{"patent_abstract":"sphaeroides"}}]}
# We  retrieved the patent number, title, and abstract.


# Config
# Change this to point to the location of the gtdbtk taxonomic classification

project_code_file <- "Project_code.txt"
MAGS_taxonomy_file <- "gtdbtk.bac120.summary.tsv"



# Required libraries

library(patentsview)
library("tidyverse")
library(stringr)
library(readr)


# Reading the project ID
# This is important to make each MAGID unique as the output for multiple projects can be the same e.g MAG??

project_id <- read.delim(file = project_code_file, header = FALSE)
project_id <- as.character(project_id$V1)

# Read MAG taxonomy file

mags <- readr::read_tsv(MAGS_taxonomy_file)

# Here we manually replaced Glutamicibacter sp004320535 with
# Arthrobacter sp. S41 to account for differences in taxonomy
# We also removed which was a spiked control added to the samples


mag_sp2 <- mags$classification %>%
  str_split(";", simplify = T) %>%
  as_tibble() %>%
  pull("V7") %>%
  str_replace_all(c("s__|_A|_C|_E|Imtechella halotolerans"), "") %>%
  str_replace("Glutamicibacter sp004320535", "Arthrobacter sp. S41") %>%
  na_if("") %>%
  as_tibble() %>%
  pull("value")

MAGID_table <- data.frame(
  "MAGID" = mags$user_genome,
  "Classification" = mags$classification,
  "Species" = mag_sp2
)

# Remove MAGs with no taxonomy classification to the species level
MAGID_table <- subset(MAGID_table, !is.na(MAGID_table$Species))

# Initilized counter
counter <- 1
while (counter <= nrow(MAGID_table)) {
  the_species <- as.character(MAGID_table[counter, 3])
  print(paste("Working on", the_species))

  Genus <- data.frame(str_split(the_species, pattern = " "))[1, ]
  Species <- data.frame(str_split(the_species, pattern = " "))[2, ]
  Genus <- as.character(Genus[1])
  Species <- as.character(Species[1])

  fields <- c("patent_number", "patent_title", "patent_abstract")

  # Creates query with specific language used in the patent system
  query_v_3 <-
    with_qfuns(
      and(
        text_all(patent_abstract = Genus),
        text_all(patent_abstract = Species)
      )
    )

  # The actual search
  tryCatch(
    expr = {
      results <- search_pv(query_v_3, all_pages = T, fields = fields)
    },
    error = function(results) {
      results <- 0
    }
  )

  results2 <- ifelse(!exists("results"), 0, results)

  if (is.list(results2)) {
    fileout_fp <- paste0(
      project_id, "_",
      as.character(MAGID_table[counter, 1]),
      "_",
      the_species, "_US_patents.txt"
    )
    
    write.table(
      x = results$data$patents,
      file = fileout_fp,
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )
    
    # readr::write_delim(
    #   x = results$data$patents,
    #   file = file_fp,
    #   delim = "\t"
    # )
  } else {
    print("No results found")
  }

  # Update counter
  counter <- counter + 1
}




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
# [1] forcats_0.5.0     stringr_1.4.1     dplyr_1.0.2       purrr_0.3.4       readr_1.3.1       tidyr_1.1.2      
# [7] tibble_3.0.3      ggplot2_3.3.5     tidyverse_1.3.0   patentsview_0.3.0
# 
# loaded via a namespace (and not attached):
# [1] Rcpp_1.0.9        cellranger_1.1.0  pillar_1.4.6      compiler_3.5.3    dbplyr_1.4.4      R.utils_2.12.0   
# [7] R.methodsS3_1.8.2 tools_3.5.3       digest_0.6.29     R.cache_0.16.0    jsonlite_1.8.0    lubridate_1.8.0  
# [13] lifecycle_1.0.1   gtable_0.3.0      pkgconfig_2.0.3   rlang_0.4.11      reprex_0.3.0      cli_3.2.0        
# [19] rstudioapi_0.13   DBI_1.1.0         curl_4.3.2        haven_2.3.1       styler_1.7.0      withr_2.5.0      
# [25] xml2_1.3.2        httr_1.4.4        fs_1.5.0          generics_0.1.3    vctrs_0.3.4       hms_0.5.3        
# [31] grid_3.5.3        tidyselect_1.1.0  glue_1.6.2        R6_2.5.1          readxl_1.3.1      rematch2_2.1.2   
# [37] modelr_0.1.8      blob_1.2.1        magrittr_2.0.3    backports_1.1.10  scales_1.1.1      ellipsis_0.3.2   
# [43] rvest_0.3.6       assertthat_0.2.1  colorspace_2.0-3  stringi_1.5.3     munsell_0.5.0     broom_0.8.0      
# [49] crayon_1.5.1      R.oo_1.25.0