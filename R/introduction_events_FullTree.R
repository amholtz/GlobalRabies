################################################################
##       Detecting Human mediated importation events full tree
##
##  author: Andrew Holtz
##  date: 21 October 2022
###############################################################

library(dplyr)
library(tidyverse)
library(ggplot2)
library(plotly)
library(treeio)
library(phangorn)
library(cepiigeodist)
library(countrycode)
library(reshape)
library(data.table)
library(DT)
library(optparse)

##Parsing


option_list = list(
  make_option(c("-m", "--meta"), type="character", default=NULL, 
              help="metadata file path", metavar="character"),
  make_option(c("-t", "--tree"), type="character", default=NULL, 
              help="pastml nexus tree file path", metavar="character"),
  make_option(c("-n", "--newick"), type="character", default=NULL, 
              help="pastml newick tree file path", metavar="character"),
  make_option(c("-a", "--annotations"), type="character", default=NULL, 
              help="pastml annotations file path", metavar="character"),
  make_option(c("-p", "--probabilities"), type="character", default=NULL, 
              help="pastml probabilities file path", metavar="character"),
  make_option(c("-x", "--prefix"), type="character", default=NULL, 
              help="prefix to identify the tree", metavar="character")
);


opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

meta <- read.delim(opt$meta)
tree <- read.nexus(opt$tree)
annotations_full <- read.delim(opt$annotations)
probabilities_full <- read.delim(opt$probabilities)


textfile <- readLines(opt$newick)


country_relations <- as.data.frame(cepiigeodist::dist_cepii)
country_relations$parent_country <- countrycode(country_relations$iso_o, 'iso3c', "country.name")
country_relations$child_country <- countrycode(country_relations$iso_d, 'iso3c', "country.name")
country_relations <- country_relations %>% select(parent_country,child_country,dist, contig)


select.tip.or.node <- function(element, tree) {
  ifelse(element < Ntip(tree)+1,
         tree$tip.label[element],
         tree$node.label[element-Ntip(tree)])
}

edge_table <- data.frame(
  "parent" = tree$edge[,1],
  "par.name" = sapply(tree$edge[,1],
                      select.tip.or.node,
                      tree = tree),
  "child" = tree$edge[,2],
  "chi.name" = sapply(tree$edge[,2],
                      select.tip.or.node,
                      tree = tree),
  "branch_length" = tree$edge.length)

edge_table$tips<- sapply(edge_table$child, function(x){
  phangorn::Descendants(tree, x, type = 'tips')
})

names(edge_table)[4] <- 'node'
names(edge_table)[2] <- 'parent_node'

#Let's see all the node character predictions and the probabilities 
prob_long_full <- melt(probabilities_full, id="node")
prob_long_full$value <- as.numeric(prob_long_full$value)
names(prob_long_full)[2] <- 'Country'
prob_long_full$Country <- gsub(" ", "", prob_long_full$Country, fixed = TRUE)
prob_long_full$Country <- gsub(".", "", prob_long_full$Country, fixed = TRUE)


#Let's join this with the annotations that have been attached for PastML
annotations_full$Country <- gsub(".", "", annotations_full$Country, fixed = TRUE)
annotations_full$Country <- gsub(" ", "", annotations_full$Country, fixed = TRUE)
prob_nodes_full <- left_join(annotations_full, prob_long_full, by = c("node", "Country"))

#Let's remove all anotations that have less than 0.5 probability 
prob_nodes_0.5_full <- prob_nodes_full %>% filter(prob_nodes_full$value > 0.5)

#join with edge_table
edge_table <- left_join(edge_table, prob_nodes_0.5_full, by = 'node')
edge_table <- arrange(edge_table,child)

names(edge_table)[4] <- 'child_node'
names(edge_table)[2] <- 'node'
names(edge_table)[7] <- 'child_country'


edge_table <- left_join(edge_table, prob_nodes_0.5_full, by = 'node')

names(edge_table)[2] <- 'parent_node'
names(edge_table)[9] <- 'parent_country'




##CREATING ANNOTATION TABLE FROM NEWICK ANNOTATED TREE
##NOTE: TO TEST VISUALIZE YOU MUST USE CAT() AND NOT SIMPLY PRINT
textfile <- gsub(pattern = '\\(', replace = "", x = textfile)
textfile  <- gsub(pattern = '\\)', replace = "", x = textfile)
textfile  <- gsub(pattern = '&&NHX:', replace = "", x = textfile)
textfile  <- gsub(pattern = '\\,', replace = "", x = textfile)
textfile  <- gsub(pattern = '\\]', replace = '\\1\n', x = textfile)
textfile  <- gsub(pattern = '\\[', replace = ':', x = textfile)
textfile  <- gsub(pattern = ':', replace = ',', x = textfile)
textfile  <- gsub(pattern = 'date=', replace = '', x = textfile)
textfile  <- gsub(pattern = 'value=', replace = '', x = textfile)
textfile  <- gsub(pattern = 'Country=', replace = '', x = textfile)

text <- as.data.frame(textfile)
fileConn<-file(paste0("../data/ACR_Results/Country/Full_Tree/", opt$prefix, "_annotations_table.csv"))
writeLines(text$textfile, fileConn)

##This file written can now be opened in R as a dataframe :)

annotations <- read_csv(paste0("../data/ACR_Results/Country/Full_Tree/", opt$prefix, "_annotations_table.csv"), col_names = FALSE)

annotations <- select(annotations, -X2, -X4)
annotations_child <- annotations
annotations_parent <- annotations
names(annotations_child)[1] <- 'child_node'
names(annotations_child)[2] <- 'child_date'

names(annotations_parent)[1] <- 'parent_node'
names(annotations_parent)[2] <- 'parent_date'

annotations <- left_join(edge_table, annotations_child, by = 'child_node')
annotations <- left_join(annotations, annotations_parent, by = 'parent_node')

annotations <- annotations %>% filter(!is.na(child_country))
annotations <- annotations %>% filter(!is.na(parent_country))

setDT(annotations)[, num_tips_child:=sapply(tips,length)]


##### DATA PREPARATION FROM TREES OVER NOW ANALYSIS

#Was there an introduction from parent to child?
annotations$intro <- ifelse(annotations$child_country == annotations$parent_country, 'False', 'True')
total_intro <- annotations %>% group_by(intro) %>% tally()
annotations <- annotations %>% filter(intro == "True")

#Distance between parent and children nodes
annotations <- left_join(annotations,country_relations, by = c("parent_country", "child_country"))
annotations$rate <- annotations$dist/annotations$branch_length
annotations <- annotations %>% filter(!is.na(dist))
#mean distance of introduction events between children and parent
print(paste0("The mean distance traveled by introduction events: ", mean(annotations$dist), ' km'))


#mean distance of introduction events between children and parent (exluding
#introductions from neighboring countries)
annotations_intro_nonneigh <- annotations %>% filter(contig == 0)
print(paste0("The number of introductions to non-neighboring countries: ", length(annotations_intro_nonneigh$child_node)))
print(paste0("The mean distance traveled over time for the introductions: ", median(annotations_intro_nonneigh$rate), ' km/year'))

annotations_intro_nonneigh <- annotations_intro_nonneigh %>% 
  select(parent_country, child_country, parent_date, child_date,
         branch_length, num_tips_child, dist, rate) %>% 
  arrange(desc(rate))

##Remove not interesting countries
not_valid <- c("Canada", "Russia", "Greenland")
print(paste0("Removing arctic countries since impossible to tell distance covered (Canada, Russia, Greenland)"))


annotations_intro_nonneigh <- annotations_intro_nonneigh %>% filter(!child_country %in% not_valid)
annotations_intro_nonneigh <- annotations_intro_nonneigh %>% filter(!parent_country %in% not_valid)

##Remove locations with distance less than 2000km (except for morocco - Europe, India-Oman, Mexico-Cuba)
annotations_intro_nonneigh_762 <- annotations_intro_nonneigh %>% filter(dist == 762.6525)
annotations_intro_nonneigh_1784 <- annotations_intro_nonneigh %>% filter(dist == 1784.223)
annotations_intro_nonneigh_1936 <- annotations_intro_nonneigh %>% filter(dist == 1936.222)
annotations_intro_nonneigh <- annotations_intro_nonneigh %>% filter(dist > 2000)
annotations_intro_nonneigh <- rbind(annotations_intro_nonneigh,
                                    annotations_intro_nonneigh_1784,
                                    annotations_intro_nonneigh_762,
                                    annotations_intro_nonneigh_1936)

annotations_intro_nonneigh <- annotations_intro_nonneigh %>% filter(rate > 200)

print(paste0("Keeping only countries with transmission over sea or ocean (Morocco-Europe, India-Oman, Mexico-Cuba), or cross distance greater than 2000km and a km travel over year is greater than 200 km/year"))



write.csv(annotations_intro_nonneigh, paste0('../data/ACR_Results/Country/Full_Tree/', opt$prefix, '_introtable.csv'), row.names = FALSE)






