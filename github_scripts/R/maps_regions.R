
################################################################
##       CREATING BASE WORLD MAP WITH COLORS DEFINED AS IN PHYLOGENETIC TREE ITOL

##  author: Andrew Holtz
##  date: 18 October 2022
###############################################################


library(cepiigeodist)
library(countrycode)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(googleVis)
library(plotly)
library(rworldmap)
library(rgdal)

setwd("/Volumes/NGS_Viroscreen/aholtz/euroME/GlobalRabies/total/alignment/better_alignment/complete_tree/complete_dating/canine")

meta <- read.delim("/Volumes/NGS_Viroscreen/aholtz/euroME/GlobalRabies/total/alignment/better_alignment/complete_tree/complete_dating/canine/rscripts/data_from_scripts/meta_16124_Reg_Exc.tab")


country_relations <- as.data.frame(cepiigeodist::dist_cepii)
countries_regions <- country_relations %>% select(iso_o) %>% unique()
countries_regions$iso_o <- str_replace(countries_regions$iso_o, 'ZAR', 'COD')
countries_regions$iso_o <- str_replace(countries_regions$iso_o, 'ROM', 'ROU')
countries_regions$region23 <- countrycode(countries_regions$iso_o, 'iso3c', 'region23')

add <- data.frame(rbind(c('SSD', 'Middle Africa'),
c('MNE', 'Southern Europe'),
c('SRB', 'Southern Europe'),
c('XXK', 'Southern Europe')))

names(add)<-c("iso_o","region23")

countries_regions <- rbind(countries_regions, add)
countries_regions <- arrange(countries_regions, region23)
countries_regions <- countries_regions %>% filter(region23 %in% meta$region23)

#Define Color Pallette

region_colors <-c(
  "Caribbean"=	"#cc8866",
  "Central America"	="#99cc99",
  "Central Asia"	="#3365cc",
  "Eastern Africa"	="#cc66aa",
  "Eastern Asia"	="#cc9933",
  "Eastern Europe"	="#66cc88",
  "Middle Africa"	="#9999cc",
  "Northern Africa"	="#cc3366",
  "Northern America"	="#ECF233",
  "Northern Europe"	="#33cc99",
  "South America"	="#8766cc",
  "South-Eastern Asia"	="#cc9999",
  "Southern Africa"	="#aacc66",
  "Southern Asia"	="#99cccc",
  "Southern Europe"	="#9933cc",
  "Western Africa"	="#66cc33",
  "Western Asia"	="#66aacc",
  "Western Europe"	="#cc99cc")



malMap <- joinCountryData2Map(countries_regions, joinCode = "ISO3",
                              nameJoinColumn = "iso_o") 
# This will join your malDF data.frame to the country map data

mapCountryData(malMap, nameColumnToPlot="region23", catMethod = "categorical",
               missingCountryCol = gray(.8), addLegend='FALSE', colourPalette = region_colors)



