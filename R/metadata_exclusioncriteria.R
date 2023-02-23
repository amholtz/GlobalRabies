####################################################
##       Defining Exclusion Criteria To Metadata at the End
##
## author: Andrew Holtz
## creation date: 05 September 2022

###################################################


library(dplyr)
library(readr)
library(treeio)
library(tidyr)
library(optparse)


`%notin%` <- Negate(`%in%`)
##Parsing

option_list = list(
  make_option(c("-d", "--meta"), type="character", default=NULL, 
              help="metadata file path", metavar="character"),
  make_option(c("-f", "--full_tree"), type="character", default=NULL, 
              help="bats/skunks/canine tree file path", metavar="character"),
  make_option(c("-c", "--canine_orig"), type="character", default=NULL, 
              help="canine tree before LSD outlier removal", metavar="character"),
  make_option(c("-l", "--canine_lsd"), type="character", default=NULL, 
              help="full canine tree post LSD outlier removal", metavar="character"),
  make_option(c("-s", "--canine_cons"), type="character", default=NULL, 
              help="consensus canine tree file path", metavar="character"),
  make_option(c("-o", "--meta_out"), type="character", default=NULL, 
              help="metadata output file path", metavar="character")
);

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

##### TEST
#opt$meta <-'/Volumes/NGS_Viroscreen/aholtz/euroME/project/GlobalRabies/data/metadata_edited.tab'
#opt$full_tree<- '/Volumes/NGS_Viroscreen/aholtz/euroME/project/GlobalRabies/data/genewise_aln_RABV.nwk'
#opt$canine_orig <-'/Volumes/NGS_Viroscreen/aholtz/euroME/project/GlobalRabies/data/RABV_canine10209.nwk'
#opt$canine_lsd <- '/Volumes/NGS_Viroscreen/aholtz/euroME/project/GlobalRabies/data/TempestRooted1327_WGSRate_OutRem.date.nwk'
#opt$canine_cons <-'/Volumes/NGS_Viroscreen/aholtz/euroME/project/GlobalRabies/data/ACR_Results/Country/Consensus_Tree/consensus_tree6096.nwk'


meta_full <- read.delim(opt$meta, na.strings = c("","NA"))
fullSpecies_tree <- read.tree(opt$full_tree)
canine_tree_orig <- read.tree(opt$canine_orig)
canine_tree_lsd <- read.tree(opt$canine_lsd)
canine_6096Tree <- read.tree(opt$canine_cons)

meta_missingdate <- meta_full %>% filter(is.na(Collection_Date))
meta_missingcountry <- meta_full %>% filter(is.na(Country))
meta_1972 <- meta_full %>% filter(Collection_Date <1972)
meta_less99 <- meta_full %>% filter(Length < 100)

meta_9663 <- Reduce(union, list(meta_missingcountry$Accession,
                                meta_missingdate$Accession,
                                meta_1972$Accession,
                                meta_less99$Accession))
meta_protcode <- meta_full %>% filter(Accession %notin% fullSpecies_tree$tip.label)
meta_notCanine <- meta_full %>% filter(Accession %notin% canine_tree_orig$tip.label)
meta_LSD2Out2 <- meta_full %>% filter(Accession %notin% canine_tree_lsd$tip.label)
meta_6096 <- meta_full %>% filter(Accession %in% canine_6096Tree$tip.label)

meta_lsd2Out <- meta_protcode %>% filter(Length > 1000)

meta_full$exclusion <- ifelse(meta_full$Accession == 'NC_001542', 'reference',
                              ifelse(meta_full$Accession %in% meta_missingdate$Accession, "no_date",
                                     ifelse(meta_full$Accession %in% meta_missingcountry$Accession, "no_country",
                                            ifelse(meta_full$Accession %in% meta_1972$Accession, "pre_1972",
                                                   ifelse(meta_full$Accession %in% meta_less99$Accession, "less_100", 
                                                          ifelse(meta_full$Accession %in% meta_lsd2Out$Accession, "LSDOutlier",
                                                                 ifelse(meta_full$Accession %in% meta_protcode$Accession, 'noncoding_region',
                                                                        ifelse(meta_full$Accession %in% meta_notCanine$Accession, 'not_canine',
                                                                               ifelse(meta_full$Accession %in% meta_LSD2Out2$Accession, 'LSDOutlier',
                                                                                      ifelse(meta_full$Accession %in% meta_6096$Accession, 'agg_included',
                                                                                             'full_included'))))))))))

meta_full$fragment <- replace_na(meta_full$fragment, "clean_omitted")
meta_full$exclusion <- ifelse(meta_full$fragment == 'removed', 'noncoding_region', meta_full$exclusion)

write_delim(meta_full, opt$meta_out, delim = '\t')
