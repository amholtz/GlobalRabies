####################################################
##       Defining Exclusion Criteria To Metadata
##
## author: Andrew Holtz
## creation date: 05 September 2022

###################################################


library(dplyr)
library(readr)
library(treeio)
library(tidyr)

setwd("/Volumes/NGS_Viroscreen/aholtz/euroME/total/alignment/better_alignment/complete_tree/complete_dating/canine")

`%notin%` <- Negate(`%in%`)


meta <- read.delim("/Volumes/NGS_Viroscreen/aholtz/euroME/total/meta_RABV_cleaned_clade_gene.tab", na.strings = c("","NA"))
meta_full <- read.delim("/Volumes/NGS_Viroscreen/aholtz/euroME/total/allRABV.tab", na.strings = c("","NA"))
fullSpecies_tree <- read.tree('/Volumes/NGS_Viroscreen/aholtz/euroME/total/alignment/better_alignment/complete_tree/complete_genewise_aln_RABV_minusINTGEN_bs0.01.nwk')
canine_tree_orig <- read.tree("/Volumes/NGS_Viroscreen/aholtz/euroME/total/alignment/better_alignment/complete_tree/complete_dating/canine/canine.nwk")
canine_tree_full <- read.tree("/Volumes/NGS_Viroscreen/aholtz/euroME/total/alignment/better_alignment/complete_tree/complete_dating/canine/dated_final_trees_Full/TempestRooted1327_WGSRate_OutRem.date.nwk")
canine_6096Tree <- read.tree("/Volumes/NGS_Viroscreen/aholtz/euroME/total/alignment/better_alignment/complete_tree/complete_dating/canine/6096_FullSubAnnotatedTrees/Final_Adapted.Fullpruned_subsampledTips_pastml/named.tree_dated.Fullpruned_subsampledTips.nwk")

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
meta_LSD2Out2 <- meta_full %>% filter(Accession %notin% canine_tree_full$tip.label)
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
meta <- meta %>% select(-Length, -Geo_Location, -Country, -Host, -Collection_Date)
meta_full <- left_join(meta_full, meta, by = 'Accession')
meta_full$fragment <- replace_na(meta_full$fragment, "clean_omitted")
meta_full$exclusion <- ifelse(meta_full$fragment == 'removed', 'noncoding_region', meta_full$exclusion)

meta_full<- meta_full %>% select(-exclusion,exclusion)

write_delim(meta_full, "/Volumes/NGS_Viroscreen/aholtz/euroME/total/meta_full_exclusion.tab", delim = '\t')

DT::datatable(meta_full, extensions = 'Buttons', options = list(
  autoWidth = TRUE, 
  pageLength = 45,
  dom = 'Bfrtip',
  buttons = c('copy', 'csv', 'excel', 'pdf', 'print')),
  filter = list(
    position = 'top', clear = FALSE))
