# Process annotation summary
# This script processes summaries tables from bakta


# This part needs to be change to match your files

project_code_file <- "Project_code.txt"
mag_counts_file <- "median_coverage_genomes.tsv"


## Required libraries
library(tidyverse)


# Reading the project ID
# This is important to make each MAGID unique as the output for multiple projects can be the same e.g MAG??

project_id <- read.delim(file = project_code_file, header = FALSE)
project_id <- as.character(project_id$V1)


## Read count file, transform the table, and add project code

count_table <- read.delim2(file=mag_counts_file, header = TRUE, dec=".")
names(count_table)[1] = "SampleID"

count_table_wide <- count_table %>%
  pivot_longer(cols = names(count_table)[2] : names(count_table)[dim(count_table)[2]],
               names_to = "MAG",
               values_to = "Median_coverage")
count_table_wide$MAGID = paste0(project_id, "_", count_table_wide$MAG) 
count_table_wide = count_table_wide[, c("MAGID", "SampleID", "Median_coverage")]

write.table(
  x = count_table_wide,
  file = "median_MAG_coverage.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
