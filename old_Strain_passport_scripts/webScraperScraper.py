#!/usr/bin/env python

#------------------------------------------------------------
# Pull webscraping data in subdirectories into database tables
#
# Environment: conda activate strainpass
# Usage: python webScraperScraper.py -d <data directory>
#------------------------------------------------------------
import os
import argparse as ap
import pandas as pd
from glob import glob
from sqlalchemy import create_engine, types


# Setup Database Engine
#------------------------------------------------------------
#Postgres
#engine = create_engine('postgresql://passport:microbial@localhost:5432/StrainPassport')
#sqlite
engine = create_engine('sqlite:///StrainPassport.db', echo=False)


# Setup Arguments
#------------------------------------------------------------
pre = "MAG" #sub-directory prefix

parser  = ap.ArgumentParser(description='scraper configuration')
parser.add_argument("-d", "--directory", type = str, help = "Path to main directory", required = True)
args    = parser.parse_args()
mainDir = getattr(args, "directory")

 
# Functions
#------------------------------------------------------------

def getStrainDirs(mainDir, pre):
    ''' Return list of strain directories'''
    strainDirs = glob(f"{mainDir}/{pre}*")
    return strainDirs

def pushStrainTable(strainDirs, tableName, pathTemplate, sep = "\t", engine = engine):
    ''' Push statistics for each strain to database table'''
    try:
        table = pd.read_sql(tableName, con=engine)
        tableExists = True
    except:
        tableExists = False

    for strainDir in strainDirs:
        strain       = os.path.split(strainDir)[1]
        strain_path  = pathTemplate.format(strainDir, strain)
        strain_stats = pd.read_csv(strain_path,sep=sep,quotechar='\'', keep_default_na=False, encoding='utf8')
        strain_stats["ID"] = int(strain.replace(pre, "")) #Table ID == MAG# with prefix removed
        strain_stats.set_index('ID', inplace = True)
        strain_stats = strain_stats.astype(str)

        #push to database
        if tableExists and table[table["ID"] == strain].empty:
            strain_stats.to_sql(tableName,con=engine,index=True,if_exists='append')
        elif not tableExists:
            strain_stats.to_sql(tableName,con=engine,index=True,if_exists='append')



# Main
#------------------------------------------------------------

def main():
    strainDirs = getStrainDirs(mainDir, pre)
    #Can interate a dictionary to reduce the repetitiveness of this section
    #Assembly Statistics
    pushStrainTable(strainDirs,
            tableName    = "AssemblyStats",
            pathTemplate = "{}/Assembly/{}.fasta.gz.assembly_stats.txt")

    #Taxonomy
    pushStrainTable(strainDirs,
            tableName    = "Taxonomy",
            pathTemplate = "{}/Taxonomy/{}_taxonomy")

    #Completeness
    pushStrainTable(strainDirs,
            tableName    = "Completeness",
            pathTemplate = "{}/Taxonomy/{}_genome_completeness")
    
    #AntibioticNCBI
    pushStrainTable(strainDirs,
            tableName    = "AntibioticNCBI",
            pathTemplate = "{}/Antibiotic/{}.fasta_ncbi.csv", sep =",")

    #AntibioticMegaRes
    pushStrainTable(strainDirs,
            tableName    = "AntibioticMegaRes",
            pathTemplate = "{}/Antibiotic/{}.fasta_megares.csv", sep =",")

    #Virulence
    pushStrainTable(strainDirs,
            tableName    = "Virulence",
            pathTemplate = "{}/Virulence/{}.fasta_vfdb.csv", sep =",")

    #Missing: Metabolism, Reference, Patents

# Main
#------------------------------------------------------------
if __name__ == '__main__':
    main()
