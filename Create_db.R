#install.packages("RSQLite")
library('RSQLite')
library(DBI)
library(tidyverse)

# Config
MAGS_taxonomy_file = "gtdbtk.bac120.summary.tsv"
project_id_file = "Project_code.txt"

mydb <- dbConnect(RSQLite::SQLite(), "my-db.sqlite")


project_id = read.delim(file=project_id_file, header=FALSE)
project_id = as.character(project_id$V1)



# List of species used for literature search

library(readr)
mags <- readr::read_tsv(MAGS_taxonomy_file)
#mags
mag_sp <- mags$classification %>%
  str_split(";", simplify = T) %>%
  as_tibble() %>%
  pull("V7") %>%
  str_replace_all(c("s__|_A|_C|_E|Imtechella halotolerans"), "") %>%
  str_replace("Glutamicibacter sp004320535", "Arthrobacter sp. S41") 
#%>%
#  na_if("") %>%
#  as_tibble() %>%
#filter(!is.na(value)) %>%
#  pull("value")
mags$mag_id = paste0(project_id, "_", mags$user_genome)
mags$Species2 = mag_sp

# Taxonomy table
tax_table = read.delim2(file="gtdbtk.bac120.summary.tsv")
tax_table$project_id = project_id
rownames(tax_table) = tax_table$user_genome 
tax_table$mag_id = paste0(tax_table$project_id, "_", tax_table$user_genome) 
class_table = matrix(nrow = length(tax_table$user_genome), ncol=8) 
rownames(class_table) = tax_table$user_genome
colnames(class_table) = c("mag_id",
                          "domain", "phylum","class","order",
                          "family", "genus","species"
                          )
class_table = as.data.frame(class_table)
class_table$mag_id = tax_table$mag_id 

for (taxa in tax_table$user_genome) {
  classification_string = as.character(tax_table[taxa, "classification"])
  classification_string = strsplit(x=classification_string, split =  ";")
  
  # Get different parts of the string
  my_domain = gsub(classification_string[[1]][1], pattern = "^d__", replacement = "")
  my_phylum = gsub(classification_string[[1]][2], pattern = "^p__", replacement = "")
  my_class = gsub(classification_string[[1]][3], pattern = "^c__", replacement = "")
  my_order = gsub(classification_string[[1]][4], pattern = "^o__", replacement = "")
  my_family = gsub(classification_string[[1]][5], pattern = "^f__", replacement = "")
  my_genus = gsub(classification_string[[1]][6], pattern = "^g__", replacement = "")
  my_species = gsub(classification_string[[1]][7], pattern = "^s__", replacement = "")
  class_table[taxa, "domain"] = my_domain
  class_table[taxa, "phylum"] = my_phylum
  class_table[taxa, "class"] = my_class
  class_table[taxa, "family"] = my_family
  class_table[taxa, "order"] = my_order
  class_table[taxa, "genus"] = my_genus
  class_table[taxa, "species"] <- my_species
  }

tax_table2 = merge(tax_table, class_table, by="mag_id")
mags = mags[, c("mag_id", "Species2")]
tax_table2 = merge(tax_table2, mags, by="mag_id")
#names(tax_table2)
taxonomy = tax_table2[, c("mag_id", "user_genome",
                         "domain", "phylum","class","order",
                          "family", "genus","species", "Species2")
                     ]

dbWriteTable(mydb, "taxonomy", taxonomy, append=T)

# Assembly

Assembly = read.delim2(file="Results/genome_summaries.tsv")
dbWriteTable(mydb, "assembly", Assembly, append=T)


# Biorxiv

Biorxiv = read.csv(file="Results/biorxiv.csv")
dbWriteTable(mydb, "Biorxiv", Biorxiv, append=T)


# Abundances

abundances = read.delim2(file="Results/median_MAG_coverage.tsv")
dbWriteTable(mydb, "abundances", abundances, append=T)

# Project_samples

project = data.frame("project_id"=project_id, "SampleID" = unique(abundances$SampleID))
dbWriteTable(mydb, "Project", project, append=T)


# Annotations

annotations = read.delim2(file="Results/merged_annotations.tsv", quote = "")
dbWriteTable(mydb, "annotations", annotations, append=T)


# Pubmed

pubmed_files = list.files(path ="Results", patter = "*all_pubmed.txt" )

#All pubmed results
pubmed_results = data.frame(matrix(ncol=13))
names(pubmed_results)= c("Species", "pmid", "doi", "title", "abstract", "year", "month", "day", "jabbrv",
                       "journal", "keywords", "lastname", "firstname")

taxonomy2 = taxonomy
taxonomy2$MAGID = taxonomy2$mag_id


counter = 1
while (counter <= length(pubmed_files)){
  pubmed = read.delim2(pubmed_files[counter], quote="")
  sp_MAGID = gsub(pubmed_files[counter], pattern = "_all_pubmed.txt", replacement = "") 
  pubmed$MAGID = sp_MAGID
  print(paste0("Woprking on ", sp_MAGID))
  pubmed2 = merge(pubmed, taxonomy2[, c("MAGID", "Species2")], by = "MAGID")
  names(pubmed2)[14] = "Species"
  pubmed2 = pubmed2[, c("Species", "pmid", "doi", "title", "abstract", "year", "month", "day", "jabbrv",
                        "journal", "keywords", "lastname", "firstname")]
  
  pubmed_results = rbind(pubmed_results, pubmed2)
  counter=counter + 1
}
pubmed_results = subset(pubmed_results, !(is.na(pmid)))

dbWriteTable(mydb, "pubmed", pubmed_results, append=T)


# Patents

patent_files = list.files(path ="Results", pattern = "*US_patents.txt" )

patent_results = data.frame(matrix(ncol=4))
names(patent_results)= c("Species", "patent_number", "patent_title", "patent_abstract")
counter = 1
while (counter <= length(patent_files)){
  patent = read.delim2(patent_files[counter], quote="")
  sp_MAGID = gsub(patent_files[counter], pattern = "_US_patents.txt", replacement = "") 
  #sp_MAGID
  patent$MAGID = sp_MAGID
  print(paste0("Woprking on ", sp_MAGID))
  patent2 = merge(patent, taxonomy2[, c("MAGID", "Species2")], by = "MAGID")
  names(patent2)[5] = "Species"
  patent2 = patent2[, c("Species", "patent_number", "patent_title", "patent_abstract")]
  patent_results = rbind(patent_results, patent2)
  counter=counter + 1
}
patent_results = subset(patent_results, !(is.na(patent_number)))

dbWriteTable(mydb, "patents", patent_results, append=T)


# ARMfinderplus

armfinder = read.csv(file="Results/armfinder_annotation.csv")
dbWriteTable(mydb, "armfinder", armfinder, append=T)

# VFDB
vfdb = read.csv(file="Results/vfdb-annotations.csv")
dbWriteTable(mydb, "vfdb", vfdb, append=T)

## BacDive
bacdive = read.csv(file="Results/bacdive_data.csv")
dbWriteTable(mydb, "Bacdive", bacdive, append=T)


RSQLite::dbDisconnect(mydb)
