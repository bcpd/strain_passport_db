
library(tidyverse)
library(readr)
library(plyr)

rm(list=ls())

library(RISmed)

api_key = "0cb41d0e0c2136aff7e0478bc4182be1cb09"


# Create list of organisms names to search

mags <- readr::read_tsv("../mixtures/gtdbtk.bac120.summary.tsv")

mag_sp <- mags$classification %>%
  str_split(";", simplify = T) %>%
  as_tibble() %>%
  pull("V7") %>%
  str_replace_all(c("s__|_A|_C|_E|Imtechella halotolerans"), "") %>%
  str_replace("Glutamicibacter sp004320535", "Arthrobacter sp. S41") %>%
  na_if("") %>% as_tibble() %>% filter(!is.na(value)) %>% pull("value")


#mag_sp

# [1] "Leuconostoc citreum"        "Bacillus inaquosorum"       "Bacillus velezensis"       
# [4] "Bacillus paralicheniformis" "Bacillus megaterium"        "Brevibacillus reuszeri"    
# [7] "Arthrobacter sp. S41"       "Pseudomonas fluorescens"    "Pseudomonas aeruginosa"    
# [10] "Cereibacter sphaeroides" 

mag_sp[7]  = "Arthrobacter S41"
#da_species = mag_sp[10]


mag_sp = c("Acinetobacter guillouiae", 
           "Acinetobacter baylyi", 
           "Arthrobacter globiformis", 
           "Bacillus subtilis", 
           "Bacillus licheniformis ", 
           "Bacillus amyloliquefaciens", 
           "Pseudomonas stutzeri", 
           "Rhodococcus rhodochrous", 
           "Rhodospirillum rubrum", 
           "Rhodobacter sphaeroides", 
           "Rhodopseudomonas palustris"
)



for (da_species in mag_sp) {

  print(paste("Working on", da_species))
  
  Genus = data.frame(str_split(da_species,  pattern = " "))[1,]
  Sp = data.frame(str_split(da_species,  pattern = " "))[2,]
  term = paste0( Genus, "[TIAB] AND ", Sp, "[TIAB]")
  #term
  
  for (query_year in seq(from = 2000, to = 2021, by =1 )){
    Sys.sleep(4)
    print(paste("Year:", query_year))
    rm(res)
  
    res <- EUtilsSummary(query = term, type="esearch", db="pubmed",
                         datetype='pdat', mindate=query_year,
                         maxdate=query_year, retmax=1000,) 
    QueryCount(res)
    
    if(QueryCount(res) > 10000){
      print("Need more records")
    }
    
    
    if(QueryCount(res) == 0) {
      print("No results found")
      
    } else {
    
      fetch <- EUtilsGet(res, type="efetch", db="pubmed")
      
      Abstracts<-data.frame('Abstract'=AbstractText(fetch))
      
      AuthorList<-Author(fetch)
      Last1 <- sapply(AuthorList, function(x) subset (x, x$order ==1))
      Last1 = (t(Last1))
      Last1 = (as.data.frame(Last1))
      FirstAU = Last1$LastName
      FirstAU <- data.frame(matrix(unlist(FirstAU), nrow=length(FirstAU), byrow=TRUE))
      
      PMIDs = data.frame('PMID'=PMID(fetch))
      Journal = Title(fetch)
      Titles = ArticleId(fetch)
      
      DOIs = DOI(fetch)
      Years = YearPubDate(fetch)
      Issue = Issue(fetch)
      Volume = Volume(fetch)
      Pages = MedlinePgn(fetch)
      
      final_results = data.frame('PubmedID' = PMIDs,
                                 'FirstAuthor' = FirstAU,
                                 "Year" = Years, "Title" = Titles, "Journal" = Journal,
                                 "Volume" = Volume, "Issue"=Issue,
                                 "Pages" = Pages, "Abstract" = Abstracts, "DOI" =DOIs)
      
      names(final_results)[2]="FirstAuthor"
      
      write_delim(x=final_results, delim = "\t", file=paste0(Genus,"_", Sp,"_", query_year, "_pubmed.txt"))
    }
  }
}




# This part is done in Linux to join files from different years

# Create folder
# mkdir tojoin

# Copy results there (change for each organism)

# cp Pseudomonas_fluorescens_20* tojoin

# Capture header of first file

# head -1 Pseudomonas_fluorescens_2000_pubmed.txt > header

# join all the files
# cat *txt >> joined

# filter the file to remove all the lines with the PMID word i.e. header lines
# grep -v 'PMID' joined > joined2

# Remove duplicated entries

# sort joined2| uniq -u > joined3

# Add header back to joined file
# cat header joined3 > Pseudomonas_fluorescens_pubmed.txt




# sessionInfo()
# R version 4.0.3 (2020-10-10)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 17763)
# 
# Matrix products: default
# 
# locale:
#   [1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252   
# [3] LC_MONETARY=English_United States.1252 LC_NUMERIC=C                          
# [5] LC_TIME=English_United States.1252    
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] plyr_1.8.6      forcats_0.5.1   stringr_1.4.0   dplyr_1.0.3     purrr_0.3.4     readr_1.4.0     tidyr_1.1.2    
# [8] tibble_3.0.5    ggplot2_3.3.3   tidyverse_1.3.0 RISmed_2.2     
# 
# loaded via a namespace (and not attached):
#   [1] Rcpp_1.0.6        cellranger_1.1.0  pillar_1.4.7      compiler_4.0.3    dbplyr_2.1.0      tools_4.0.3      
# [7] jsonlite_1.7.2    lubridate_1.7.9.2 lifecycle_0.2.0   gtable_0.3.0      pkgconfig_2.0.3   rlang_0.4.10     
# [13] reprex_1.0.0      cli_2.2.0         DBI_1.1.1         rstudioapi_0.13   haven_2.3.1       withr_2.4.1      
# [19] xml2_1.3.2        httr_1.4.2        fs_1.5.0          generics_0.1.0    vctrs_0.3.6       hms_1.0.0        
# [25] grid_4.0.3        tidyselect_1.1.0  glue_1.4.2        R6_2.5.0          fansi_0.4.2       readxl_1.3.1     
# [31] modelr_0.1.8      magrittr_2.0.1    backports_1.2.1   scales_1.1.1      ellipsis_0.3.1    rvest_0.3.6      
# [37] assertthat_0.2.1  colorspace_2.0-0  stringi_1.5.3     munsell_0.5.0     broom_0.7.4       crayon_1.3.4 