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
library(cepiigeodist)
library(countrycode)
library(ggplot2)
library(seqinr)
library(readr)
library(dplyr)
library(stringr)
library(data.table)
library(taxize)
library(optparse)

##Parsing

option_list = list(
  make_option(c("-d", "--meta"), type="character", default=NULL, 
              help="metadata file path", metavar="character"),
  make_option(c("-t", "--host_table"), type="character", default=NULL, 
              help="host species table path", metavar="character"),
  make_option(c("-a", "--aln"), type="character", default=NULL, 
              help="alignment file path", metavar="character"),
  make_option(c("-n", "--out_n_text"), type="character", default=NULL, 
              help="output N gene file path", metavar="character"),
  make_option(c("-p", "--out_p_text"), type="character", default=NULL, 
              help="output P gene file path", metavar="character"),
  make_option(c("-m", "--out_m_text"), type="character", default=NULL, 
              help="output M gene file path", metavar="character"),
  make_option(c("-g", "--out_g_text"), type="character", default=NULL, 
              help="output G gene file path", metavar="character"),
  make_option(c("-l", "--out_l_text"), type="character", default=NULL, 
              help="output L gene file path", metavar="character"),
  make_option(c("-w", "--out_wgs_text"), type="character", default=NULL, 
              help="output wgs gene file path", metavar="character"),
  make_option(c("-o", "--meta_out"), type="character", default=NULL, 
              help="metadata output file path", metavar="character")
);


opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)


##
meta <- read.delim(opt$meta)

## Creating a simpler column for clade definitions
meta$clade_simple <- ifelse(str_detect(meta$Clade, "Bats"), "Bat_Clade", 
                            ifelse(str_detect(meta$Clade, "Africa_2"), "Africa2_Clade", 
                                   ifelse(str_detect(meta$Clade, "Africa_3"), "Africa3_Clade", 
                                          ifelse(str_detect(meta$Clade, "Arctic"), "Arctic_Clade", 
                                                 ifelse(str_detect(meta$Clade, "Cosmo"), "Cosmopolitan_Clade", 
                                                        ifelse(str_detect(meta$Clade, "Asian"), "Asian_Clade", meta$Clade))))))


## Creating a simpler column for host species definitions 
host_table <- read_csv(opt$host_table)
meta$host_simple <- host_table$host_simple[match(meta$Host, host_table$Host)]

## Creating a column for regional definitions

meta$region23<- countryname(meta$Country, destination = "region23")

## Creating a column for colonial history

country_relations <- as.data.frame(cepiigeodist::dist_cepii)
country_relations$parent_country <- countrycode(country_relations$iso_o, 'iso3c', "country.name")
country_relations$child_country <- countrycode(country_relations$iso_d, 'iso3c', "country.name")

#Define colonial powers:
colonial_powers <- c("United Kingdom", "France", "Spain", "Portugal", "Russia")

country_relations <- country_relations %>% select(parent_country,child_country,colony)
country_relations$parent_country <- str_replace(country_relations$parent_country, 'United States', 'USA')
country_relations$child_country <- str_replace(country_relations$child_country, 'United States', 'USA')

country_relations <- country_relations %>% filter(parent_country %in% colonial_powers) %>% 
  filter(colony == 1)

UK_colony <- country_relations %>% filter(parent_country == "United Kingdom") %>% select(child_country) %>%
  rbind(c("United Kingdom")) %>% 
  rbind(c("China"))

FR_colony <- country_relations %>% filter(parent_country == "France") %>% select(child_country) %>% rbind(c("France"))
RU_colony <- country_relations %>% filter(parent_country == "Russia") %>% select(child_country) %>% rbind(c("Russia"))
PR_colony <- country_relations %>% filter(parent_country == "Portugal") %>% select(child_country) %>% rbind(c("Portugal"))
SP_colony <- country_relations %>% filter(parent_country == "Spain") %>% select(child_country) %>% rbind(c("Spain"))

meta$colony <- ifelse(meta$Country %in% UK_colony$child_country, 'British Empire', 
                      ifelse(meta$Country %in% FR_colony$child_country, 'French Empire',
                             ifelse(meta$Country %in% RU_colony$child_country, 'Russian Empire',
                                    ifelse(meta$Country %in% SP_colony$child_country, 'Spainish Empire',
                                           ifelse(meta$Country %in% PR_colony$child_country, 'Portugese Empire', 'Non-colonial')))))


## Creating text files for sequence accession numbers that contain nucleotides in a subgenomic region

aln = read.fasta(opt$aln)

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



#######
meta$fragment <- ifelse(meta$Length > 10000, 'WGS',
                        ifelse(meta$Accession %in% Ngenes, 'N',
                               ifelse(meta$Accession %in% Pgenes, 'P',
                                      ifelse(meta$Accession %in% Mgenes, 'M',
                                             ifelse(meta$Accession %in% Ggenes, 'G',
                                                    ifelse(meta$Accession %in% Lgenes, 'L', 'removed'))))))

write_delim(meta, opt$meta_out, delim = '\t', quote = 'none')

## Removing sequences with data quality issues
meta <- meta %>% filter(!is.na(Country)) %>% filter(!is.na(Collection_Date)) %>% 
  filter(Collection_Date > 1971) %>% filter(Length > 99)

meta_n <- meta %>% filter(meta$fragment == 'N') %>% select(Accession) %>% 
  write.csv(opt$out_n_text,
            row.names = FALSE, quote = FALSE, col.names = NA)
meta_p <- meta %>% filter(meta$fragment == 'P') %>% select(Accession) %>% 
  write.csv(opt$out_p_text,
            row.names = FALSE, quote = FALSE, col.names = FALSE)
meta_m <- meta %>% filter(meta$fragment == 'M') %>% select(Accession) %>% 
  write.csv(opt$out_m_text,
            row.names = FALSE, quote = FALSE, col.names = FALSE)
meta_g <- meta %>% filter(meta$fragment == 'G') %>% select(Accession) %>% 
  write.csv(opt$out_g_text,
            row.names = FALSE, quote = FALSE, col.names = FALSE)
meta_l <- meta %>% filter(meta$fragment == 'L') %>% select(Accession) %>% 
  write.csv(opt$out_l_text,
            row.names = FALSE, quote = FALSE, col.names = FALSE)
meta_wgs <- meta %>% filter(meta$fragment == 'WGS') %>% select(Accession) %>% 
  write.csv(opt$out_wgs_text,
            row.names = FALSE, quote = FALSE, col.names = FALSE)



