################################################################
##       Colonizations PASTML Preparation & PhyCova Matrix 
##
##  author: Andrew Holtz
##  date: 30 August 2022
###############################################################
## The plan is to make a column in data that shows whether the sequence
## Is coming from a country that used to be part of a major colonial power 
## Between 1600 and 1950. The colonial powers here are the largest 6 during that time:
## British Empire
## Russian Empire + NOT USSR SATEILLTE STATES
## Qing dynasty
## Spanish Empire
## French Empire
## Portugese Empire


library(cepiigeodist)
library(countrycode)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(plotly)

meta <- read.delim("/Volumes/NGS_Viroscreen/aholtz/euroME/total/alignment/better_alignment/complete_tree/complete_dating/canine/rscripts/data_from_scripts/meta_16124_Reg_Exc.tab")

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
write_delim(meta, "/Volumes/NGS_Viroscreen/aholtz/euroME/total/meta_RABV_cleaned_clade_gene.tab", delim = '\t')

## PHYCOVA Matrix Prep



meta_included <- meta %>% filter(exclusion == "agg_included")

country <- meta_included %>% group_by(Country) %>% tally()
colony <- meta_included %>% group_by(colony) %>% tally()

colony$colony <- factor(colony$colony, levels = 
                          colony$colony[order(-colony$n)])

colony_colors <- c("British Empire"=	"#7fc97f",
                  "French Empire"	="#beaed4",
                  "Non-colonial"	="#fdc086",
                  "Portugese Empire"	="#ffff99",
                  "Russian Empire"	="#386cb0",
                  "Spainish Empire"	="#f0027f")

colony <- ggplot(colony, aes(colony, n, fill = colony)) + geom_bar(stat = "identity") +
  labs(title = "RABV Sequences by Colonizing Power", x = "Colonizing Power", y =  "Count") +
  scale_fill_manual(values = colony_colors)+
  theme_bw() +
  theme_classic() 

country <- select(country, Country, colony)

country$colony <- ifelse(country$Country %in% UK_colony$child_country, 'British Empire', 
                         ifelse(country$Country %in% FR_colony$child_country, 'French Empire',
                                ifelse(country$Country %in% RU_colony$child_country, 'Russian Empire',
                                       ifelse(country$Country %in% SP_colony$child_country, 'Spainish Empire',
                                              ifelse(country$Country %in% PR_colony$child_country, 'Portugese Empire', 'Non-colonial')))))


matrix <- read.csv("~/Dropbox/B.1.214/phylo/glm_data/distances.csv", header=FALSE)
matrix$V1 <- gsub(" ", "", matrix$V1, fixed = TRUE)
matrix$V2 <- gsub(" ", "", matrix$V2, fixed = TRUE)
matrix <- matrix %>% select(V1, V2) %>% filter(V1 %in% country$Country) %>% filter(V2 %in% country$Country)
matrix<- left_join(matrix, country, by = c("V1" = "Country"))
matrix<- left_join(matrix, country, by = c("V2" = "Country"))
names(matrix)[1] <- 'country1'
names(matrix)[2] <- 'country2'
names(matrix)[3] <- 'colony1'
names(matrix)[4] <- 'colony2'
matrix$same_colony <- ifelse(matrix$colony1 == matrix$colony2, matrix$colony2, 0)
matrix$same_colony <- ifelse(matrix$same_colony == 'Non-colonial', 0, matrix$same_colony)
matrix$same_colony <- ifelse(matrix$country1 == matrix$country2, 1, matrix$same_colony)
matrix$same_colony <- ifelse(matrix$same_colony >= 1, 1, 0)

un1 <- unlist(unique(matrix[1:2]))
matrix[1:2] <- lapply(matrix[1:2], factor, levels = un1)
matrix_glm <- as.data.frame.matrix(xtabs(matrix$same_colony ~ matrix$country1 + matrix$country2, matrix))

write.csv(matrix_glm, '/Volumes/NGS_Viroscreen/aholtz/euroME/total/alignment/better_alignment/complete_tree/complete_dating/canine/rscripts/data_from_scripts/country_colony_matrix.csv', quote = FALSE)





