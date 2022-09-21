
#####
library(BacDive)
#modified from: https://github.com/jlab/algorithm_traits_microbiome/blob/a9fa8c12b34bff8a20581e5e6568cea1944577a9/bacdive_custom_func.R

# Using BacDive requires free registration

bacdive <- open_bacdive("padimitriu@gmail.com", "phase2_spdb")

#get BacDive IDs
request_bac <- function(taxon){
  taxon_id <- request(object = bacdive, # request returns just the bacdive IDs 
                      query = taxon,     
                      search = "taxon", page= 0L, handler = NULL, sleep = 0.1)
  message(ifelse(taxon_id$count == 1,taxon_id$results, 
                 paste0(taxon," has ", taxon_id$count, " strains!")))
  return(unlist(taxon_id$results))
}

# Obtain species from gtdb table
#example
# taxa <- c("Pseudomonas aeruginosa", "Bacillus cereus", "Arthrobacter globiformis")

taxa_list <- map(taxa, request_bac)


## Retrieve Bacdive data


bacdat <- function(ID){

x <- fetch(bacdive, ID)
y <- x$results
z <- unlist(y)
 
# for antibiotic resistance data  
dfpa <- data.frame(map_dfr(names(z), as.data.frame),
           map_dfr(unname(z), as.data.frame))
dfpa <- dfpa %>%
  rename(catnames = 1,
         catvalues = 2)

dat.merged <- dfpa %>%
  dplyr::group_by(catnames) %>%
  dplyr::summarise(catvalues = paste(catvalues, collapse = ";"))

arc <- dat.merged %>%
  filter(catnames == paste0(ID,".Physiology and metabolism.antibiotic resistance.ChEBI")) %>%
  select(catvalues) %>% pull()

arm <- dat.merged %>%
  filter(catnames == paste0(ID,".Physiology and metabolism.antibiotic resistance.metabolite")) %>%
  select(catvalues) %>% pull()


  df <- data.frame(
    Genus = ifelse(paste0(ID,".Name and taxonomic classification.genus")%in% names(z),z[[paste0(ID,".Name and taxonomic classification.genus")]],"unclassified"),
    
    Species = ifelse(paste0(ID,".Name and taxonomic classification.species")%in% names(z),z[[paste0(ID,".Name and taxonomic classification.species")]],"unclassified"),
    
    Synonyms = ifelse(paste0(ID,".Name and taxonomic classification.LPSN.synonyms.synonym")%in% names(z),z[[paste0(ID,".Name and taxonomic classification.LPSN.synonyms.synonym")]],NA),
    
    ncbi_taxid = ifelse(paste0(ID,".General.NCBI tax id.NCBI tax id")%in% names(z),z[[paste0(ID,".General.NCBI tax id.NCBI tax id")]], NA),
    
    taxid_matching_level = ifelse(paste0(ID,".General.NCBI tax id.Matching level")%in% names(z),z[[paste0(ID,".General.NCBI tax id.Matching level")]],NA),
    
    type_strain = ifelse(paste0(ID,".Name and taxonomic classification.type strain")%in% names(z),z[[paste0(ID,".Name and taxonomic classification.type strain")]],NA),
    
    medium_name = ifelse(paste0(ID,".Culture and growth conditions.culture medium.name")%in% names(z),z[[paste0(ID,".Culture and growth conditions.culture medium.name")]],NA),
    
    medium_link = ifelse(paste0(624,".Culture and growth conditions.culture medium.link")%in% names(z),z[[paste0(624,".Culture and growth conditions.culture medium.link")]],NA),
    
    strain_history = ifelse(paste0(ID,".General.strain history2")%in% names(z),z[[paste0(ID,".General.strain history2")]],NA),
    
    gram = ifelse(paste0(ID,".Morphology.cell morphology.gram stain")%in% names(z),z[[paste0(ID,".Morphology.cell morphology.gram stain")]],NA),
    
    length = ifelse(paste0(ID,".Morphology.cell morphology.cell length")%in% names(z),z[[paste0(ID,".Morphology.cell morphology.cell length")]],NA),
    
    width = ifelse(paste0(ID,".Morphology.cell morphology.cell width")%in% names(z),z[[paste0(ID,".Morphology.cell morphology.cell width")]],NA),
    
    motility = ifelse(paste0(ID,".Morphology.cell morphology.motility")%in% names(z),z[[paste0(ID,".Morphology.cell morphology.motility")]],NA),
    
    oxygen_tolerance = ifelse(paste0(ID,".Physiology and metabolism.oxygen tolerance.oxygen tolerance")%in% names(z),z[[paste0(ID,".Physiology and metabolism.oxygen tolerance.oxygen tolerance")]],NA),
    
    spore = ifelse(paste0(ID,".Physiology and metabolism.spore formation.spore formation")%in% names(z),z[[paste0(ID,".Physiology and metabolism.spore formation.spore formation")]],NA),
    
    aggregation_score = ifelse(paste0(ID,".Physiology and metabolism.observation.observation")%in% names(z),z[[paste0(ID,".Physiology and metabolism.observation.observation")]],NA),
    
    GC_content = ifelse(paste0(ID,".Sequence information.GC content.GC-content")%in% names(z),z[[paste0(ID,".Sequence information.GC content.GC-content")]],NA),
    
    antibiotic_resistance_chebi = ifelse(paste0(ID,".Physiology and metabolism.antibiotic resistance.ChEBI")%in% names(z),arc,NA),
    
    antibiotic_resistance_metabolite = ifelse(paste0(ID,".Physiology and metabolism.antibiotic resistance.metabolite")%in% names(z),arm,NA),  
    
    row.names = NULL,stringsAsFactors = FALSE,check.names = T)
  
  return(df)
}


bacdive_data <- map_dfr(unlist(taxa_list), bacdat) %>%
  distinct(ncbi_taxid, .keep_all = T)

write_csv(bacdive_data, "bacdive_data.csv")
