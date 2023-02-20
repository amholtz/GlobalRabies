###################################################
##       MetaData -> Number of countries by sequence choice
##
## author: Andrew Holtz
## creation date: 10 October 2022

###################################################


library(dplyr)
library(readr)
library(treeio)
library(tidyr)
library(ggplot2)
library(stringr)
library(plyr)
library(countrycode)

setwd("/Volumes/NGS_Viroscreen/aholtz/euroME/GlobalRabies/total/alignment/better_alignment/complete_tree/complete_dating/canine")

meta <- read.delim("/Volumes/NGS_Viroscreen/aholtz/euroME/GlobalRabies/total/alignment/better_alignment/complete_tree/complete_dating/canine/rscripts/data_from_scripts/meta_16124_Reg_Exc.tab")

meta_full <- meta %>% filter(exclusion == "full_included" | exclusion == "agg_included")
meta_full <- meta_full %>% group_by(fragment, Country) %>% tally()

meta_g <- meta_full %>% filter(str_detect(fragment, "G")) %>% ungroup() %>%  select(Country) %>% unique()
meta_n <- meta_full %>% filter(str_detect(fragment, "N")) %>% ungroup() %>%  select(Country) %>% unique()
meta_p <- meta_full %>% filter(str_detect(fragment, "P")) %>% ungroup() %>%  select(Country) %>% unique()
meta_l <- meta_full %>% filter(str_detect(fragment, "L")) %>% ungroup() %>%  select(Country) %>% unique()
meta_m <- meta_full %>% filter(str_detect(fragment, "M")) %>% ungroup() %>%  select(Country) %>% unique()
meta_wgs <- meta_full %>% filter(str_detect(fragment, "NPMGL")) %>% ungroup() %>%  select(Country) %>% unique()


meta_full <- meta_full %>% ungroup() %>% select(Country) %>% unique()
