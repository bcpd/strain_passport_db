# Process annotation summary
# This script processes summaries tables from bakta


# This part needs to be change to match your files

project_code_file = "Project_code.txt"
MAGS_bakta_folder = "MAGS_bakta"


## Required libraries
library(tidyverse)



# Reading the project ID
# This is important to make each MAGID unique as the output for multiple projects can be the same e.g MAG??

project_id = read.delim(file=project_code_file, header=FALSE)
project_id = as.character(project_id$V1)





# Parsing list of individual MAG annotation summaries from bakta

baktas2 = list.files(path = MAGS_bakta_folder, pattern = "MAG[0-9]+.txt", full.names = T )

bakta_stats = as.data.frame(matrix(ncol=3))
names(bakta_stats) = c("Variable", "Value", "MAGID")

for (bakta_file in baktas2){
  
  myfile = read_file(bakta_file)
  myfile2 = strsplit(x =  myfile, split = "\n")
  myfile2 = as.data.frame(myfile2)
  names(myfile2)[1] = "Variable"
  myfile2 = subset(myfile2, myfile2$Variable != "")
  myfile3 = separate(myfile2, col = Variable, sep = ":", into=c("Variable", "Value"))
  myfile3 = subset(myfile3, myfile3$Value != "")
  MAG_name = gsub(bakta_file, pattern = "[a-zA-Z0-9_/]+/", replacement="")
  MAG_name = gsub(MAG_name, pattern = ".txt$", replacement="")
  myfile3$MAGID = paste0(project_id, "_", MAG_name )
  bakta_stats =rbind(bakta_stats, myfile3)
}

# This line subsets the previous table, ignoring the first row which has NA for
# all its values
bakta_stats = bakta_stats [2:nrow(bakta_stats), c("MAGID", "Variable", "Value")]  

bakta_stats_wide = pivot_wider(data = bakta_stats, names_from = Variable, values_from = Value)
write.table(x = bakta_stats_wide, file ="genome_summaries.tsv", quote = FALSE, row.names = FALSE,
            sep = "\t" )


## Session information

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
# [1] Rcpp_1.0.9       pillar_1.4.6     compiler_3.5.3   cellranger_1.1.0 dbplyr_1.4.4     tools_3.5.3     
# [7] jsonlite_1.8.0   lubridate_1.8.0  lifecycle_1.0.1  gtable_0.3.0     pkgconfig_2.0.3  rlang_0.4.11    
# [13] reprex_0.3.0     cli_3.2.0        rstudioapi_0.13  DBI_1.1.0        haven_2.3.1      withr_2.5.0     
# [19] xml2_1.3.2       httr_1.4.4       fs_1.5.0         generics_0.1.3   vctrs_0.3.4      hms_0.5.3       
# [25] grid_3.5.3       tidyselect_1.1.0 glue_1.6.2       R6_2.5.1         readxl_1.3.1     modelr_0.1.8    
# [31] blob_1.2.1       magrittr_2.0.3   backports_1.1.10 scales_1.1.1     ellipsis_0.3.2   rvest_0.3.6     
# [37] assertthat_0.2.1 colorspace_2.0-3 stringi_1.5.3    munsell_0.5.0    broom_0.8.0      crayon_1.5.1    
