

r = getOption("repos")
r["CRAN"] = "http://cran.us.r-project.org"
options(repos = r)
#install.packages('prettydoc')
#install.packages('kableExtra')
#install.packages('knitr')
#install.packages('circlize')
#install.packages('plotly')
#install.packages('easyPubMed')
#install.packages('pathview')
# install.packages('rmdformats')

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("rtracklayer")
BiocManager::install("pathview")
#install.packages("rtracklayer")

#source("https://bioconductor.org/biocLite.R")
#biocLite("rtracklayer")

#Libraries
#---------------------------------------------------------------
library("DBI")
library("stringr")
library("rmarkdown")
library("ggplot2")
library("pander")
library("prettydoc")
library("kableExtra")
library("knitr")
library("rtracklayer")
library("circlize")
library("plotly")
library("easyPubMed")
library("pathview")
library("jsonlite")
library("tidyverse")

'%!in%' <- function(x,y)!('%in%'(x,y))

#Directories
#---------------------------------------------------------------
mainDir = paste0(getwd(), "/", "Results") #location of files necessary for generating reports
reportDir = "StrainPassports"
setwd(mainDir)

#Create Directory of Individual Library Results
dir.create(file.path(mainDir, reportDir),
           showWarnings = FALSE)


#Database
#---------------------------------------------------------------
# Create an ephemeral in-memory RSQLite database
# con <- dbConnect(RSQLite::SQLite(), "StrainPassport.db")
# tableList <- dbListTables(con)


#Data
#---------------------------------------------------------------
#Load Annotations
annotations <- read_tsv("merged_annotations.tsv")
                         
annotations$pfam_hits <- paste0(data.frame(str_split_fixed(annotations$pfam_hits, ";",2))[,1], "...")
annotations$ko_id <-   data.frame(str_split_fixed(annotations$ko_id, ",",2))[,1]

annotations <- annotations %>%
  inner_join(tibble(value = mag_sp, MAGID = mag_list), by = "MAGID")
  
#head(annotations)

#specialized KEGG modules 
# download and save json files in mainDir
#https://www.genome.jp/kegg-bin/get_htext?ko01504
#https://www.genome.jp/kegg-bin/get_htext?htext=ko02042&orgs=hsa%20eco
#https://www.genome.jp/kegg-bin/get_htext?ko02022+K07770

ko01504.resistome <- read_json("ko01504.json")
ko02042.toxins <- read_json("ko02042.json")
ko02022.twoComp <- read_json("ko02022.json")

ko01504.resistome.DF <- ko01504.resistome %>%
  as_tibble() %>%
  unnest_auto(children) %>%
  unnest_auto(children) %>%
  unnest_auto(children) %>%
  unnest_auto(children) %>%
  unnest_auto(children) %>%
  unnest_auto(children) %>%
  unnest_auto(children)

ko02042.toxins.DF <- ko02042.toxins %>%
  as_tibble() %>%
  unnest_auto(children) %>%
  unnest_auto(children) %>%
  unnest_auto(children) %>%
  unnest_auto(children) %>%
  unnest_auto(children) %>%
  unnest_auto(children) %>%
  unnest_auto(children)

ko02022.twoComp.DF <- ko02022.twoComp %>%
  as_tibble() %>%
  unnest_auto(children) %>%
  unnest_auto(children) %>%
  unnest_auto(children) %>%
  unnest_auto(children) %>%
  unnest_auto(children) 

ko02042.toxins.DF_t3 <- ko02042.toxins.DF %>%
  filter(`name...2` %in% c("Type III toxins: Intracellular toxins"))

ko02042.toxins.DF_t3$Description <- data.frame(str_split_fixed(ko02042.toxins.DF_t3$children, "  ", 2))[,2]
ko02042.toxins.DF_t3$KO          <- data.frame(str_split_fixed(ko02042.toxins.DF_t3$children, "  ", 2))[,1]
ko02042.toxins.DF_t3$children <- NULL

ko02042.toxins.DF_t12 <- ko02042.toxins.DF %>%
  filter(`name...2` %in% c("Type II toxins: Membrane damaging toxins",
                           "Type I toxins: Toxins that act from the cell surface",
                           "Toxins that damage the extracellular matrix",
                           "Not specified"))
ko02042.toxins.DF_t12$Description <- data.frame(str_split_fixed(ko02042.toxins.DF_t12$name...4, "  ", 2))[,2]
ko02042.toxins.DF_t12$KO          <- data.frame(str_split_fixed(ko02042.toxins.DF_t12$name...4, "  ", 2))[,1]
ko02042.toxins.DF_t12$children <- NULL

# colbind toxins
ko02042.toxins.DF <- bind_rows(ko02042.toxins.DF_t12, ko02042.toxins.DF_t3)


ko01504.resistome.DF$KO <- gsub("  .*", "", ko01504.resistome.DF$name...5) #Antimicrobial Resistance Genes
ko02042.toxins.DF$KO <- gsub("  .*", "", ko02042.toxins.DF$name...4) #Toxin Genes
ko02022.twoComp.DF$KO <- gsub("  .*", "", ko02022.twoComp.DF$name...4) #Two Component System Genes



#Summarize KEGG hits per genome
#use this: 
#https://www.genome.jp/kegg-bin/get_htext?ko00001.keg ; download json file

ko00001 <- read_json("ko00001.json")
KO1 <- ko00001 %>%
  as_tibble() %>%
  unnest_auto(children) %>%
  unnest_auto(children) %>%
  unnest_auto(children) %>% 
  unnest_auto(children) %>%
  unnest_auto(children) %>%
  unnest_auto(children) %>%
  unnest_auto(children) %>%
  unnest_auto(children) %>%
  unnest_auto(children) %>%
  unnest_auto(children) 

write_delim(KO1, "KO_Orthology_ko00001.txt", delim = "\t")

KO <- read_delim("KO_Orthology_ko00001.txt",
                 delim ="\t")

KO$Description <- data.frame(str_split_fixed(KO$children, "  ", 2))[,2]
KO$KO          <- data.frame(str_split_fixed(KO$children, "  ", 2))[,1]

KO.loci <- KO[match(annotations$ko_id, KO$KO), ]
#head(KO.loci)
KO.loci$Strain = annotations$MAGID


KO.level1 <- aggregate(KO.loci, by=list(KO.loci$Strain, KO.loci$name...2), FUN=length)
KO.level1 <- KO.level1[,1:3]
colnames(KO.level1) <- c("Strain", "Category", "N.genes")
#unique(KO.level1$Category)
KO.level1 <- KO.level1 %>%
  inner_join(tibble(value = mag_sp, Strain = mag_list), by = "Strain")


KO.level2 <- aggregate(KO.loci, by=list(KO.loci$Strain, KO.loci$name...3),
                       FUN=length)
KO.level2 <- KO.level2[,1:3]
colnames(KO.level2) <- c("Strain", "Category", "N.genes")
KO.level2 <- KO.level2 %>%
  inner_join(tibble(value = mag_sp, Strain = mag_list), by = "Strain")

#as.character(unique(KO.level2$Category))

# virulence genes
KO.vir <- KO.loci[which(KO.loci$name...3 == "09171 Infectious disease: bacterial"), ]
virulence.annotations <- annotations[which(annotations$ko_id %in% KO.vir$KO), ]


level2.catagories <- c("09101 Carbohydrate metabolism", "09102 Energy metabolism", "09103 Lipid metabolism",
  "09104 Nucleotide metabolism", "09105 Amino acid metabolism" , "09107 Glycan biosynthesis and metabolism",
  "09111 Xenobiotics biodegradation and metabolism" , "09121 Transcription", "09122 Translation",
  "09124 Replication and repair", "09131 Membrane transport","09132 Signal transduction",
  "09141 Transport and catabolism", "09142 Cell motility", "09143 Cell growth and death",
  "09159 Environmental adaptation", "09171 Infectious disease: bacterial", "09175 Drug resistance: antimicrobial")



#gff file
#gff <- data.frame(import.gff("annotation/genes.gff"))
#gff$Strain <- data.frame(str_split_fixed(gff$seqnames, "_", 3))[,1]
#head(gff)



#Strain list
strains <- readr::read_tsv("gtdbtk.bac120.summary.txt",
                           show_col_types = F)


# Species-level MAGs
sps <- strains$classification %>%
  str_split(";", simplify = T) %>%
  as_tibble() %>%
  pull("V7") %>%
  str_replace_all(c("s__|_A|_C|_E|Imtechella halotolerans"), "") %>%
  str_replace("Glutamicibacter sp004320535", "Arthrobacter sp. S41") %>%
  na_if("") %>% as_tibble() %>%
  bind_cols(strains) %>%
  filter(!is.na(value))

mag_list <- sps %>% pull(user_genome)
mag_sp <- sps %>%
  pull(value)
names(mag_sp) <- mag_sp



# assembly stats
genome_summaries <- read_tsv("genome_summaries.tsv") %>%
  mutate(MAGID = str_replace(MAGID, "TEST0001_", "")) %>%
  # subset to species-level MAGS
  filter(MAGID %in% mag_list)

# reorder shit
ridx <- match(mag_list, genome_summaries$MAGID)
genome_summaries <- genome_summaries[ridx,]

genome_summaries <- genome_summaries %>%
  bind_cols(mag_sp) %>%
  rename(value = 27)

genome_summaries <- genome_summaries[c("value", "MAGID", "Length", "Count", "GC", "N50", "coding density", "CDSs")]
colnames(genome_summaries) <- c("value", "MAGID", "Size (bp)", "N.contigs", "GC (%)", "N50 (bp)", "Coding Density", "N.CDSs")

# genome completeness
genome_completeness <- read_tsv("genome_completeness.tsv") %>%
  # subset to species-level MAGS
  filter(`Bin Id` %in% mag_list)

ridx2 <- match(mag_list, genome_completeness$`Bin Id`)
genome_completeness <- genome_completeness[ridx2,]

genome_completeness <- genome_completeness %>%
  bind_cols(mag_sp) %>%
  rename(value = 15)

# amr
amr <- read_csv("armfinder_annotation.csv")  %>%
  #mutate(MAGID = str_replace(MAGID, "TEST0001_", "")) %>%
  inner_join(tibble(tibble(value = mag_sp, SRA = mag_list)), by = "SRA") %>%
  rename(class = class.y)

# vfdb
vfdb <- read_csv("vfdb-annotations.csv") 
vfdbcts <- read_csv("vfdb_map.csv")

vfdb <- vfdb %>%
  rename(VFCID = VFC) %>%
  inner_join(vfdbcts, by = "VFCID") %>%
  distinct(locus, .keep_all = T)
vfdb <- vfdb %>%
  select(SRA, VFG.x, VFCID, VFcategory, VFSubcategory, psc.uniref90_id) %>%
  inner_join(tibble(value = mag_sp, SRA = mag_list), by = "SRA")

#BacDive table

bdtable <- read_csv("bacdive_data.csv") %>%
  group_by(Species) %>%
  #count() %>%
  slice_sample()


#Create pages from RMD template
#---------------------------------------------------------------
#Pass libraries to template file
#testsp = mag_sp[c(1,2)]

for(strain in mag_sp){ #[c(2,3,15)] for testing
  rmarkdown::render("StrainPassportTemplatev2.Rmd",
         output_file=file.path(mainDir, reportDir, paste0(strain,"_Passport.html")),
         params=list(new_title=paste("Strain Passport -", strain)))
}


################################################################
### Finalize and end
#dbDisconnect(con)
