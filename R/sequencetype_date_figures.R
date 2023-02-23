###################################################
##       Creation of Figures showing Sequencing by gene fragment over time
##
## author: Andrew Holtz
## creation date: 2022/09/15

###################################################


library(dplyr)
library(readr)
library(treeio)
library(tidyr)
library(ggplot2)
library(scales)
library(countrycode)

setwd("/Volumes/NGS_Viroscreen/aholtz/euroME/total/alignment/better_alignment/complete_tree/complete_dating/canine")


meta <- read.delim("/Volumes/NGS_Viroscreen/aholtz/euroME/total/alignment/better_alignment/complete_tree/complete_dating/canine/rscripts/data_from_scripts/meta_16124_Reg_Exc.tab")
canine_tree <- read.nexus("/Volumes/NGS_Viroscreen/aholtz/euroME/total/alignment/better_alignment/complete_tree/complete_dating/canine/WGSDating_CanineFull_PastML/named.tree_TempestRooted1327_WGSRate_OutRem.date.nexus")

meta_date_gene <- meta %>% group_by(Collection_Date, fragment) %>% tally()
meta$year <- as.integer(substr(meta$Collection_Date, 1, 4))

meta$continent <- countrycode(sourcevar = meta[, "Country"],
                              origin = "country.name",
                              destination = "continent")
meta <- meta %>% filter(fragment != 'removed')

fragment_colors <- c('G' = "#a6cee3", "GL" =	'#1f78b4', "L" = '#b2df8a', "M" = '#33a02c', "N" = '#fb9a99', "NPMG" = '#e31a1c', "P" = "#ff7f00", "NPMGL" = '#fdbf6f')				
host_colors <-c("Bovidae" =	"#cc7b66",
                "Canidae"	="#b299cc",
                "Felidae"	="#ccbd66",
                "Herpestidae"	="#cc99c4",
                "Hominidae"	="#c1cc33",
                "Mephitidae"	="#cc66a7",
                "Mustelidae"	="#66a0cc",
                "dog"	="#6d66cc",
                "fox"	="#66cc74",
                "other" = "#B8B8B8")

## Reorder fragments from largest to smallest
meta$fragment <- factor(meta$fragment, 
                               levels = c('GL',
                                          'NPMG',
                                          'L',
                                          'M',
                                          'P',
                                          'NPMGL',
                                          'G', 
                                          'N'))


meta_canine <- meta %>% filter(Accession %in% canine_tree$tip.label)

ggplot(data = meta) +
  geom_histogram(mapping = aes(year, fill = fragment), binwidth = 2) +
  labs(x = "Date", y = "Count", fill = "Fragment") +
  scale_fill_manual(values = fragment_colors)+
  theme_bw() +
  theme_classic() +
  ylim(0, 3200)

ggplot(data = meta) +
  geom_freqpoly(mapping = aes(year, colour = fragment), size = 2, binwidth = 2) +
  labs(x = "Date", y = "Count", fill = "Fragment") +
  scale_colour_manual(values = fragment_colors)+
  theme_bw() +
  theme_classic()

meta_continent <- meta %>% filter(fragment != 'removed') %>% group_by(continent, fragment) %>% tally() 

ggplot(data = meta_continent) +
  geom_bar(mapping = aes(x=continent, y= n, fill = fragment), position = 'dodge', stat = 'identity') +
  labs(x = "Date", y = "Count", fill = "Fragment") +
  scale_fill_manual(values = fragment_colors)+
  ylim(0, 3200) +
  theme_bw() +
  theme_classic() 

#Now everything together in a facet wrap

ggplot(data = meta) +
  geom_histogram(mapping = aes(year, fill = fragment), binwidth = 2) +
  labs(x = "Date of extraction", y = "Number of sequences", fill = "Fragment") +
  facet_wrap(~continent) +
  scale_fill_manual(values = fragment_colors)+
  theme_bw() +
  theme_classic() 

#Canine Tree

ggplot(data = meta_canine) +
  geom_histogram(mapping = aes(year, fill = fragment), binwidth = 2) +
  labs(x = "Date of extraction", y = "Number of sequences", fill = "Fragment") +
  facet_wrap(~continent) +
  scale_fill_manual(values = fragment_colors)+
  theme_bw() +
  theme_classic() +
  ylim(0, 1050)

#Facet wrap for species and continent 

top_host <- c("dog",
              "Canidae",
              "fox",
              "Bovidae",
              "Hominidae",
              "Mustelidae",
              "Felidae",
              'Mephitidae',
              "Herpestidae")

meta_canine$top_host <- ifelse(meta_canine$host_simple %in% top_host, meta_canine$host_simple, 'other')

## Reorder top host from largest to smallest
meta_canine$top_host <- factor(meta_canine$top_host, 
                                    levels = c('Herpestidae',
                                               'Felidae', 
                                               'Mephitidae',
                                               'Mustelidae',
                                               'Hominidae',
                                               'other',
                                               'Bovidae',
                                               'fox',
                                               'Canidae', 
                                               'dog'))

ggplot(data = meta_canine) +
  geom_histogram(mapping = aes(year, fill = top_host), binwidth = 2) +
  labs(x = "Date of extraction", y = "Number of sequences", fill = "Host") +
  #facet_wrap(~continent) +
  scale_fill_manual(values = host_colors)+
  theme_bw() +
  theme_classic()

## Now continents separetly

meta_canine_asia <- meta_canine %>% filter(continent == 'Asia')
meta_canine_africa <- meta_canine %>% filter(continent == 'Africa')
meta_canine_europe <- meta_canine %>% filter(continent == 'Europe')
meta_canine_america <- meta_canine %>% filter(continent == 'Americas')


#Asia - HOST
ggplot(data = meta_canine_asia) +
  geom_histogram(mapping = aes(year, fill = top_host), binwidth = 2) +
  labs(x = "Date of extraction", y = "Number of sequences", fill = "Host") +
  scale_fill_manual(values = host_colors)+
  theme_bw() +
  facet_wrap(~continent) +
  theme_classic() +
  ylim(0, 1050)

ggsave(
  "/Volumes/NGS_Viroscreen/aholtz/euroME/total/alignment/better_alignment/complete_tree/complete_dating/canine/rscripts/data_from_scripts/asia_host.pdf",
  device = 'pdf')


#Americas - HOST
ggplot(data = meta_canine_america) +
  geom_histogram(mapping = aes(year, fill = top_host), binwidth = 2) +
  labs(x = "Date of extraction", y = "Number of sequences", fill = "Host") +
  scale_fill_manual(values = host_colors)+
  theme_bw() +
  facet_wrap(~continent) +
  theme_classic() +
  ylim(0, 400)

ggsave(
  "/Volumes/NGS_Viroscreen/aholtz/euroME/total/alignment/better_alignment/complete_tree/complete_dating/canine/rscripts/data_from_scripts/americas_host.pdf",
  device = 'pdf')

#Europe - HOST
ggplot(data = meta_canine_europe) +
  geom_histogram(mapping = aes(year, fill = top_host), binwidth = 2) +
  labs(x = "Date of extraction", y = "Number of sequences", fill = "Host") +
  scale_fill_manual(values = host_colors)+
  theme_bw() +
  facet_wrap(~continent) +
  theme_classic() +
  ylim(0, 400)

ggsave(
  "/Volumes/NGS_Viroscreen/aholtz/euroME/total/alignment/better_alignment/complete_tree/complete_dating/canine/rscripts/data_from_scripts/europe_host.pdf",
  device = 'pdf')

#Africa - HOST
ggplot(data = meta_canine_africa) +
  geom_histogram(mapping = aes(year, fill = top_host), binwidth = 2) +
  labs(x = "Date of extraction", y = "Number of sequences", fill = "Host") +
  scale_fill_manual(values = host_colors)+
  theme_bw() +
  facet_wrap(~continent) +
  theme_classic() +
  ylim(0, 400)

ggsave(
  "/Volumes/NGS_Viroscreen/aholtz/euroME/total/alignment/better_alignment/complete_tree/complete_dating/canine/rscripts/data_from_scripts/africa_host.pdf",
  device = 'pdf')

#Asia - FRAGMENT
ggplot(data = meta_canine_asia) +
  geom_histogram(mapping = aes(year, fill = fragment), binwidth = 2) +
  labs(x = "Date of extraction", y = "Number of sequences", fill = "Fragment") +
  facet_wrap(~continent) +
  scale_fill_manual(values = fragment_colors)+
  theme_bw() +
  theme_classic() +
  ylim(0, 1050)

ggsave(
  "/Volumes/NGS_Viroscreen/aholtz/euroME/total/alignment/better_alignment/complete_tree/complete_dating/canine/rscripts/data_from_scripts/asia_fragment.pdf",
  device = 'pdf')

#America - FRAGMENT
ggplot(data = meta_canine_america) +
  geom_histogram(mapping = aes(year, fill = fragment), binwidth = 2) +
  labs(x = "Date of extraction", y = "Number of sequences", fill = "Fragment") +
  facet_wrap(~continent) +
  scale_fill_manual(values = fragment_colors)+
  theme_bw() +
  theme_classic() +
  ylim(0, 400)

ggsave(
  "/Volumes/NGS_Viroscreen/aholtz/euroME/total/alignment/better_alignment/complete_tree/complete_dating/canine/rscripts/data_from_scripts/america_fragment.pdf",
  device = 'pdf')

#Europe - FRAGMENT
ggplot(data = meta_canine_europe) +
  geom_histogram(mapping = aes(year, fill = fragment), binwidth = 2) +
  labs(x = "Date of extraction", y = "Number of sequences", fill = "Fragment") +
  facet_wrap(~continent) +
  scale_fill_manual(values = fragment_colors)+
  theme_bw() +
  theme_classic() +
  ylim(0, 400)

ggsave(
  "/Volumes/NGS_Viroscreen/aholtz/euroME/total/alignment/better_alignment/complete_tree/complete_dating/canine/rscripts/data_from_scripts/europe_fragment.pdf",
  device = 'pdf')

#Africa - FRAGMENT
ggplot(data = meta_canine_africa) +
  geom_histogram(mapping = aes(year, fill = fragment), binwidth = 2) +
  labs(x = "Date of extraction", y = "Number of sequences", fill = "Fragment") +
  facet_wrap(~continent) +
  scale_fill_manual(values = fragment_colors)+
  theme_bw() +
  theme_classic() +
  ylim(0, 400)

ggsave(
  "/Volumes/NGS_Viroscreen/aholtz/euroME/total/alignment/better_alignment/complete_tree/complete_dating/canine/rscripts/data_from_scripts/africa_fragment.pdf",
  device = 'pdf')

