####################################################
##                           CLEAN ALL RABV METADATA
## ####   Removing Sequences with no country, date data, and small fragments
## ####   Defining Simple Clade information
## ####   Defining host definitions
## ####   Defining the genetic region of each sequence
## ####   
## author: Andrew Holtz
## creation date: 2022/03/16
###################################################

rm(list=ls())
library(dplyr)
library(ggplot2)
library(seqinr)
library(readr)
library(dplyr)
library(stringr)
library(data.table)
library(taxize)


setwd("/Volumes/NGS_Viroscreen/aholtz/euroME")
source("/Volumes/NGS_Viroscreen/aholtz/euroME/R/helperFunctions.R")


##
meta <- read_csv("total/allRABV.csv")
clade_metadata <- read.delim("/Volumes/NGS_Viroscreen/aholtz/euroME/total/clade_metadata.txt")

clade_metadata <- select(clade_metadata, sequence.sequenceID, alignment.name)
names(clade_metadata)[1] <- 'Accession'
names(clade_metadata)[2] <- 'Clade'

meta <- meta %>% filter(!is.na(Country)) %>% filter(!is.na(Collection_Date)) %>% 
  filter(Collection_Date > 1971) %>% filter(Length > 99)

meta <- left_join(meta, clade_metadata, by = 'Accession')

meta$clade_simple <- ifelse(str_detect(meta$Clade, "Bats"), "Bat_Clade", 
                            ifelse(str_detect(meta$Clade, "Africa_2"), "Africa2_Clade", 
                                   ifelse(str_detect(meta$Clade, "Africa_3"), "Africa3_Clade", 
                                          ifelse(str_detect(meta$Clade, "Arctic"), "Arctic_Clade", 
                                                 ifelse(str_detect(meta$Clade, "Cosmo"), "Cosmopolitan_Clade", 
                                                        ifelse(str_detect(meta$Clade, "Asian"), "Asian_Clade", meta$Clade))))))

meta$species_simple <- ifelse(new$host == 'Canis lupus familiaris', "dog",
                              ifelse(new$host == 'Vulpes vulpes', "fox",
                                     ifelse(new$host == 'Nyctereutes procyonoides', "raccoon dog",
                                            ifelse(new$host == 'Canidae', new$species.y, 'other'))))

host_species <- unique(meta$Host)

x <- taxize::tax_name(host_species, get = 'family', db = 'ncbi')
family <- x
names(family)[2] <- 'species'
family <- family[-1]
family$order <- taxize::tax_name(family$family, get = 'order', db = 'ncbi')

family$keep <- ifelse(family$species == 'Canis lupus familiaris', 'dog',
                      ifelse(family$species == 'Canis lupus', 'wolf',
                             ifelse(family$species == 'Vulpes vulpes', 'fox',
                                    ifelse(family$order$order == 'Chiroptera', 'bat',
                                           family$family))))

meta$host_simple <- ifelse(meta$Host %in% family$species, family$keep, "other")


##

#aln = read.fasta("/Volumes/NGS_Viroscreen/aholtz/euroME/total/with_keeplength_RABV.fasta")
aln = read.fasta("/Users/aholtz/Downloads/with_keeplength_RABV.fasta")

# Convert alignment
aln = as.alignment(
  nb = length(aln), 
  nam = sapply(aln, function(x) attributes(x)$name),
  seq = lapply(aln, function(x) {
    attributes(x) = NULL
    x
  })
)

# Trim alignment to positions each gene position (or WGS)
aln$seqWG = lapply(aln$seq, function(x) x[72:min(length(x), 11000)])
aln$seqN = lapply(aln$seq, function(x) x[72:min(length(x), 1424)])
aln$seqP = lapply(aln$seq, function(x) x[1512:min(length(x), 2405)])
aln$seqM = lapply(aln$seq, function(x) x[2494:min(length(x), 3102)])
aln$seqG = lapply(aln$seq, function(x) x[3315:min(length(x), 4889)])
aln$seqL = lapply(aln$seq, function(x) x[5408:min(length(x), 11000)])




# Verify that all sequences contain the N gene
Ngenes = unlist(sapply(aln$nam, function(x) { 
  s = aln$seqN[[x]]
  if (sum(s %in% "-") < 1300) {
    gsub("_.*$", "", x)
  } else {
    NULL
  }
}
))

# Verify that all sequences contain the P gene
Pgenes = unlist(sapply(aln$nam, function(x) { 
  s = aln$seqP[[x]]
  if (sum(s %in% "-") < 890) {
    gsub("_.*$", "", x)
  } else {
    NULL
  }
}
))

# Verify that all sequences contain the M gene
Mgenes = unlist(sapply(aln$nam, function(x) { 
  s = aln$seqM[[x]]
  if (sum(s %in% "-") < 605) {
    gsub("_.*$", "", x)
  } else {
    NULL
  }
}
))

# Verify that all sequences contain the G gene
Ggenes = unlist(sapply(aln$nam, function(x) { 
  s = aln$seqG[[x]]
  if (sum(s %in% "-") < 1570) {
    gsub("_.*$", "", x)
  } else {
    NULL
  }
}
))

# Verify that all sequences contain the L gene
Lgenes = unlist(sapply(aln$nam, function(x) { 
  s = aln$seqL[[x]]
  if (sum(s %in% "-") < 5590) {
    gsub("_.*$", "", x)
  } else {
    NULL
  }
}
))

meta$fragment <- ifelse(meta$Length > 10000, 'WGS',
                        ifelse(meta$Accession %in% Ngenes, 'NGene',
                               ifelse(meta$Accession %in% Pgenes, 'PGene',
                                      ifelse(meta$Accession %in% Mgenes, 'MGene',
                                             ifelse(meta$Accession %in% Ggenes, 'GGene',
                                                    ifelse(meta$Accession %in% Lgenes, 'LGene', 'removed'))))))

write.csv(meta, "/Volumes/NGS_Viroscreen/aholtz/euroME/total/meta_RABV_cleaned_clade_gene.txt",
          row.names = FALSE, quote = FALSE, col.names = FALSE)
write_delim(meta, "/Volumes/NGS_Viroscreen/aholtz/euroME/total/meta_RABV_cleaned_clade_gene.tab", delim = '\t', quote = 'none')

#date file format
meta %>% select(Accession, Collection_Date) %>% 
  write.csv("/Volumes/NGS_Viroscreen/aholtz/euroME/total/meta_RABV_cleaned_clade_gene_dates_treetime.txt",
            row.names = FALSE, quote = FALSE, col.names = FALSE)

meta_n <- meta %>% filter(meta$fragment == 'NGene') %>% select(Accession) %>% 
  write.csv("/Volumes/NGS_Viroscreen/aholtz/euroME/total/meta_n.txt",
            row.names = FALSE, quote = FALSE, col.names = FALSE)
meta_p <- meta %>% filter(meta$fragment == 'PGene') %>% select(Accession) %>% 
  write.csv("/Volumes/NGS_Viroscreen/aholtz/euroME/total/meta_p.txt",
            row.names = FALSE, quote = FALSE, col.names = FALSE)
meta_m <- meta %>% filter(meta$fragment == 'MGene') %>% select(Accession) %>% 
  write.csv("/Volumes/NGS_Viroscreen/aholtz/euroME/total/meta_m.txt",
            row.names = FALSE, quote = FALSE, col.names = FALSE)
meta_g <- meta %>% filter(meta$fragment == 'GGene') %>% select(Accession) %>% 
  write.csv("/Volumes/NGS_Viroscreen/aholtz/euroME/total/meta_g.txt",
            row.names = FALSE, quote = FALSE, col.names = FALSE)
meta_l <- meta %>% filter(meta$fragment == 'LGene') %>% select(Accession) %>% 
  write.csv("/Volumes/NGS_Viroscreen/aholtz/euroME/total/meta_l.txt",
            row.names = FALSE, quote = FALSE, col.names = FALSE)
meta_wgs <- meta %>% filter(meta$fragment == 'WGS') %>% select(Accession) %>% 
  write.csv("/Volumes/NGS_Viroscreen/aholtz/euroME/total/meta_wgs.txt",
            row.names = FALSE, quote = FALSE, col.names = FALSE)










