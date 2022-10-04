# Patent search


# We searched for patents using the patentsview package in R.
# This packages allow for searchs in the US patent system using its
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

MAGS_taxonomy_file = "gtdbtk.bac120.summary.tsv"



# Required libraries

library(patentsview)
library("tidyverse")
library(stringr)
library(readr)


mags <- readr::read_tsv(MAGS_taxonomy_file)

# Here we manually replaced Glutamicibacter sp004320535 with
# Arthrobacter sp. S41 to account for differences in taxonomy

mag_sp <- mags$classification %>%
  str_split(";", simplify = T) %>%
  as_tibble() %>%
  pull("V7") %>%
  str_replace_all(c("s__|_A|_C|_E|Imtechella halotolerans"), "") %>%
  str_replace("Glutamicibacter sp004320535", "Arthrobacter sp. S41") %>%
  na_if("") %>%
  as_tibble() %>%
  filter(!is.na(value)) %>%
  pull("value")

# This is the expected value for the list of species (mag_sp)
#
# >mag_sp
# [1] "Leuconostoc citreum"        "Bacillus inaquosorum"       "Bacillus velezensis"
# [4] "Bacillus paralicheniformis" "Bacillus megaterium"        "Brevibacillus reuszeri"
# [7] "Arthrobacter sp. S41"       "Pseudomonas fluorescens"    "Pseudomonas aeruginosa"
# [10] "Cereibacter sphaeroides"



for (the_species in mag_sp) {
  print(paste("Working on", the_species))

  Genus <- data.frame(str_split(the_species, pattern = " "))[1, ]
  Species <- data.frame(str_split(the_species, pattern = " "))[2, ]

  fields <- c("patent_number", "patent_title", "patent_abstract")

  # Creates query with specific language used in the patent system
  query_v_3 <-
    with_qfuns(
      and(
        text_all(patent_abstract = Genus),
        text_all(patent_abstract = Species)
      )
    )

  # rm(results)
  # rm(results2)

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
    readr::write_delim(
      x = results$data$patents,
      file = paste0(the_species, "_US_patents.txt"), delim = "\t"
    )
  } else {
    print("No results found")
  }
}


# sessionInfo()
# R version 4.0.3 (2020-10-10)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 17763)
#
# locale:
#   [1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252
#   [3] LC_MONETARY=English_United States.1252 LC_NUMERIC=C
#   [5] LC_TIME=English_United States.1252
#
#   attached base packages:
#     [1] stats     graphics  grDevices utils     datasets  methods   base
#
#   other attached packages:
#     [1] stringr_1.4.0     dplyr_1.0.3       plyr_1.8.6        tidyr_1.1.2       patentsview_0.2.2 RISmed_2.2
#
#   loaded via a namespace (and not attached):
#     [1] Rcpp_1.0.6       rstudioapi_0.13  magrittr_2.0.1   hms_1.0.0        tidyselect_1.1.0 R6_2.5.0
#   [7] rlang_0.4.10     fansi_0.4.2      tools_4.0.3      cli_2.2.0        DBI_1.1.1        ellipsis_0.3.1
#   [13] assertthat_0.2.1 tibble_3.0.5     lifecycle_0.2.0  crayon_1.3.4     purrr_0.3.4      readr_1.4.0
#   [19] vctrs_0.3.6      glue_1.4.2       stringi_1.5.3    compiler_4.0.3   pillar_1.4.7     generics_0.1.0
#   [25] pkgconfig_2.0.3
