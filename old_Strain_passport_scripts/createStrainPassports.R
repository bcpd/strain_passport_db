#Create strain passports using data from sqlite databae

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

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("rtracklayer")
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

'%!in%' <- function(x,y)!('%in%'(x,y))

#Directories
#---------------------------------------------------------------
mainDir = "C:/Users/dlevyboo/Documents/MBI/StrainPassport"
reportDir = "StrainPassports"
setwd(mainDir)

#Create Directory of Individual Library Results
dir.create(file.path(mainDir, reportDir), 
           showWarnings = FALSE)


#Database
#---------------------------------------------------------------
# Create an ephemeral in-memory RSQLite database
con <- dbConnect(RSQLite::SQLite(), "StrainPassport.db")
tableList <- dbListTables(con)


#Data
#---------------------------------------------------------------
#Load Annotations
annotations <- read.delim("annotation/annotations.tsv",
                          stringsAsFactors = FALSE,
                          sep = "\t",
                          header = TRUE)
annotations$pfam_hits <- paste0(data.frame(str_split_fixed(annotations$pfam_hits, ";",2))[,1], "...")
annotations$kegg_id <-   data.frame(str_split_fixed(annotations$kegg_id, ",",2))[,1]
head(annotations)

#specialized KEGG modules
ko01504.resistome.DF    <- read.table("annotation/ko01504.resistome.tsv", sep="\t", header = T)
ko02042.toxins.DF       <- read.table("annotation/ko02042.toxins.tsv", sep="\t", header = T)
ko02022.twoComp.DF      <- read.table("annotation/ko02022.twoComponent.tsv", sep="\t", header = T)

ko01504.resistome.DF$KO <- gsub("  .*", "", ko01504.resistome.DF$name4) #Antimicrobial Resistance Genes
ko02042.toxins.DF$KO <- gsub("  .*", "", ko02042.toxins.DF$name4) #Toxin Genes
ko02022.twoComp.DF$KO <- gsub("  .*", "", ko02022.twoComp.DF$name3) #Two Component System Genes



#Summerize KEGG hits per genome
KO <- read.delim("annotation/KO_Orthology_ko00001.txt",
                 header = T, sep="\t")
KO$Description <- data.frame(str_split_fixed(KO$KO, "  ", 2))[,2]
KO$KO          <- data.frame(str_split_fixed(KO$KO, "  ", 2))[,1]

KO.loci <- KO[match(annotations$kegg_id, KO$KO), ]
head(KO.loci)
KO.loci$Strain = annotations$fasta

KO.level1 <- aggregate(KO.loci, by=list(KO.loci$Strain, KO.loci$Level1), FUN=length)
KO.level1 <- KO.level1[,1:3]
colnames(KO.level1) <- c("Strain", "Catagory", "N.genes")
unique(KO.level1$Catagory)

KO.level2 <- aggregate(KO.loci, by=list(KO.loci$Strain, KO.loci$Level2), FUN=length)
KO.level2 <- KO.level2[,1:3]
colnames(KO.level2) <- c("Strain", "Catagory", "N.genes")
as.character(unique(KO.level2$Catagory))

KO.vir <- KO.loci[which(KO.loci$Level2 == "09171 Infectious disease: bacterial"), ]
virulence.annotations <- annotations[which(annotations$kegg_id %in% KO.vir$KO), ]

level2.catagories <- c("09101 Carbohydrate metabolism", "09102 Energy metabolism", "09103 Lipid metabolism", 
  "09104 Nucleotide metabolism", "09105 Amino acid metabolism" , "09107 Glycan biosynthesis and metabolism", 
  "09111 Xenobiotics biodegradation and metabolism" , "09121 Transcription", "09122 Translation",
  "09124 Replication and repair", "09131 Membrane transport","09132 Signal transduction", 
  "09141 Transport and catabolism", "09142 Cell motility", "09143 Cell growth and death",                            
  "09159 Environmental adaptation", "09171 Infectious disease: bacterial", "09175 Drug resistance: antimicrobial")       



#gff file
gff <- data.frame(import.gff("annotation/genes.gff"))
gff$Strain <- data.frame(str_split_fixed(gff$seqnames, "_", 3))[,1]
head(gff)





#Strain list
strains <- read.table("strains.txt", stringsAsFactors = F)$V1


#Create pages from RMD template
#---------------------------------------------------------------
#Pass libraries to template file
for(strain in strains){ #[c(2,3,15)] for testing
  render("StrainPassportTemplate.Rmd",
         output_file=file.path(mainDir, reportDir, paste0(strain,"_Passport.html")),
         params=list(new_title=paste("Strain Passport -", strain)))
}


################################################################
### Finalize and end
dbDisconnect(con)


