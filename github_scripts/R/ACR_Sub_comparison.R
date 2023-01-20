################################################################
##       Comparison of Nodes in ACR of Canine RABV Tree
##
##  # Compare internal nodes in subsampled trees to full-tree
##  # Review PastML ACR (MPPA) on country-level for each resolved node between subsampled and full-tree
##  # Create new consensus tree with union of sequences and keep ACR state if shared between the subsamples
##  # Create new annotation file for PastML HTML Compressed Tree Visualization
##
##  author: Andrew Holtz
##  date: 25 July 2022
###############################################################

library(dplyr)
library(tidyverse)
library(ggplot2)
library(plotly)
library(treeio)
library(phangorn)fgdh
library(cepiigeodist)
library(countrycode)
library(reshape)
library(data.table)
library(DT)

setwd("/Volumes/NGS_Viroscreen/aholtz/euroME/GlobalRabies/total/alignment/better_alignment/complete_tree/complete_dating/canine")

meta <- read.delim("/Volumes/NGS_Viroscreen/aholtz/euroME/GlobalRabies/total/meta_RABV_cleaned_clade_gene.tab")
tree_full <- read.nexus("/Volumes/NGS_Viroscreen/aholtz/euroME/GlobalRabies/total/alignment/better_alignment/complete_tree/complete_dating/canine/WGSDating_CanineFull_PastML/named.tree_TempestRooted1327_WGSRate_OutRem.date.nexus")
annotations_full <- read.delim("/Volumes/NGS_Viroscreen/aholtz/euroME/GlobalRabies/total/alignment/better_alignment/complete_tree/complete_dating/canine/WGSDating_CanineFull_PastML/combined_ancestral_states.tab")
probabilities_full <- read.delim("/Volumes/NGS_Viroscreen/aholtz/euroME/GlobalRabies/total/alignment/better_alignment/complete_tree/complete_dating/canine/WGSDating_CanineFull_PastML/marginal_probabilities.character_Country.model_F81.tab")


tree_sub1 <- read.nexus("/Volumes/NGS_Viroscreen/aholtz/euroME/GlobalRabies/total/alignment/better_alignment/complete_tree/complete_dating/canine/smart_sub/finished_subsampled/WGSDating_5000_1_PastML/named.tree_TempEstRooted_subsampled_5000_1_OutRem_WGSRate_SavedFromNexus.nexus")
annotations_sub1 <- read.delim("/Volumes/NGS_Viroscreen/aholtz/euroME/GlobalRabies/total/alignment/better_alignment/complete_tree/complete_dating/canine/smart_sub/finished_subsampled/WGSDating_5000_1_PastML/combined_ancestral_states.tab")
probabilities_sub1 <- read.delim("/Volumes/NGS_Viroscreen/aholtz/euroME/GlobalRabies/total/alignment/better_alignment/complete_tree/complete_dating/canine/smart_sub/finished_subsampled/WGSDating_5000_1_PastML/marginal_probabilities.character_Country.model_F81.tab")

tree_sub2 <- read.nexus("/Volumes/NGS_Viroscreen/aholtz/euroME/GlobalRabies/total/alignment/better_alignment/complete_tree/complete_dating/canine/smart_sub/finished_subsampled/WGSDating_5000_2_PastML/named.tree_Subsample5000_2_TempestRooted1327_WGSRate_OutRem.date.nexus")
annotations_sub2 <- read.delim("/Volumes/NGS_Viroscreen/aholtz/euroME/GlobalRabies/total/alignment/better_alignment/complete_tree/complete_dating/canine/smart_sub/finished_subsampled/WGSDating_5000_2_PastML/combined_ancestral_states.tab")
probabilities_sub2 <- read.delim("/Volumes/NGS_Viroscreen/aholtz/euroME/GlobalRabies/total/alignment/better_alignment/complete_tree/complete_dating/canine/smart_sub/finished_subsampled/WGSDating_5000_2_PastML/marginal_probabilities.character_Country.model_F81.tab")


tree_sub3 <- read.nexus("/Volumes/NGS_Viroscreen/aholtz/euroME/GlobalRabies/total/alignment/better_alignment/complete_tree/complete_dating/canine/smart_sub/finished_subsampled/WGSDating_5000_3_PastML/named.tree_TempEstRooted_subsampled_5000_3_OutRem_SavedFromNexus_LSD2.nexus")
annotations_sub3 <- read.delim("/Volumes/NGS_Viroscreen/aholtz/euroME/GlobalRabies/total/alignment/better_alignment/complete_tree/complete_dating/canine/smart_sub/finished_subsampled/WGSDating_5000_3_PastML/combined_ancestral_states.tab")
probabilities_sub3 <- read.delim("/Volumes/NGS_Viroscreen/aholtz/euroME/GlobalRabies/total/alignment/better_alignment/complete_tree/complete_dating/canine/smart_sub/finished_subsampled/WGSDating_5000_3_PastML/marginal_probabilities.character_Country.model_F81.tab")


tree_sub4 <- read.nexus("/Volumes/NGS_Viroscreen/aholtz/euroME/GlobalRabies/total/alignment/better_alignment/complete_tree/complete_dating/canine/smart_sub/finished_subsampled/WGSDating_5000_4_PastML/named.tree_TempEstRooted_subsampled_5000_4_OutRem_SavedFromNexus_LSD2.nexus")
annotations_sub4 <- read.delim("/Volumes/NGS_Viroscreen/aholtz/euroME/GlobalRabies/total/alignment/better_alignment/complete_tree/complete_dating/canine/smart_sub/finished_subsampled/WGSDating_5000_4_PastML/combined_ancestral_states.tab")
probabilities_sub4 <- read.delim("/Volumes/NGS_Viroscreen/aholtz/euroME/GlobalRabies/total/alignment/better_alignment/complete_tree/complete_dating/canine/smart_sub/finished_subsampled/WGSDating_5000_4_PastML/marginal_probabilities.character_Country.model_F81.tab")


tree_sub5 <- read.nexus("/Volumes/NGS_Viroscreen/aholtz/euroME/GlobalRabies/total/alignment/better_alignment/complete_tree/complete_dating/canine/smart_sub/finished_subsampled/WGSDating_5000_5_PastML/named.tree_TempEstRooted_subsampled_5000_5_OutRem_SavedFromNexus_LSD2.nexus")
annotations_sub5 <- read.delim("/Volumes/NGS_Viroscreen/aholtz/euroME/GlobalRabies/total/alignment/better_alignment/complete_tree/complete_dating/canine/smart_sub/finished_subsampled/WGSDating_5000_5_PastML/combined_ancestral_states.tab")
probabilities_sub5 <- read.delim("/Volumes/NGS_Viroscreen/aholtz/euroME/GlobalRabies/total/alignment/better_alignment/complete_tree/complete_dating/canine/smart_sub/finished_subsampled/WGSDating_5000_5_PastML/marginal_probabilities.character_Country.model_F81.tab")

####### Remove tips in subsamples that are not in final tree

'%!in%' <- function(x,y)!('%in%'(x,y))

sub1_ni_full <- tree_sub1$tip.label %>% as.data.frame()
sub1_ni_full$notin <- tree_sub1$tip.label %!in% tree_full$tip.label
sub1_ni_full <- sub1_ni_full %>% filter(notin == TRUE)
tree_sub1 <- drop.tip(tree_sub1, sub1_ni_full$.)

sub2_ni_full <- tree_sub2$tip.label %>% as.data.frame()
sub2_ni_full$notin <- tree_sub2$tip.label %!in% tree_full$tip.label
sub2_ni_full <- sub2_ni_full %>% filter(notin == TRUE)
tree_sub2 <- drop.tip(tree_sub2, sub2_ni_full$.)

sub3_ni_full <- tree_sub3$tip.label %>% as.data.frame()
sub3_ni_full$notin <- tree_sub3$tip.label %!in% tree_full$tip.label
sub3_ni_full <- sub3_ni_full %>% filter(notin == TRUE)
tree_sub3 <- drop.tip(tree_sub3, sub3_ni_full$.)

sub4_ni_full <- tree_sub4$tip.label %>% as.data.frame()
sub4_ni_full$notin <- tree_sub4$tip.label %!in% tree_full$tip.label
sub4_ni_full <- sub4_ni_full %>% filter(notin == TRUE)
tree_sub4 <- drop.tip(tree_sub4, sub4_ni_full$.)

sub5_ni_full <- tree_sub5$tip.label %>% as.data.frame()
sub5_ni_full$notin <- tree_sub5$tip.label %!in% tree_full$tip.label
sub5_ni_full <- sub5_ni_full %>% filter(notin == TRUE)
tree_sub5 <- drop.tip(tree_sub5, sub5_ni_full$.)

#######


country_relations <- as.data.frame(cepiigeodist::dist_cepii)
country_relations$parent_country <- countrycode(country_relations$iso_o, 'iso3c', "country.name")
country_relations$child_country <- countrycode(country_relations$iso_d, 'iso3c', "country.name")
country_relations <- country_relations %>% select(parent_country,child_country,dist, contig)


####
select.tip.or.node <- function(element, tree) {
  ifelse(element < Ntip(tree)+1,
         tree$tip.label[element],
         tree$node.label[element-Ntip(tree)])
}


edge_table_full <- data.frame(
  "child" = tree_full$edge[,2],
  "node" = sapply(tree_full$edge[,2],
                  select.tip.or.node,
                  tree = tree_full))

#Let's see all the node character predictions and the probabilities 
prob_long_full <- melt(probabilities_full, id="node")
prob_long_full$value <- as.numeric(prob_long_full$value)
names(prob_long_full)[2] <- 'Country'
meta <- read.delim("/Volumes/NGS_Viroscreen/aholtz/euroME/GlobalRabies/total/meta_RABV_cleaned_clade_gene.tab")
prob_long_full$Country <- gsub(" ", "", prob_long_full$Country, fixed = TRUE)
prob_long_full$Country <- gsub(".", "", prob_long_full$Country, fixed = TRUE)


#Let's join this with the annotations that have been attached for PastML
annotations_full$Country <- gsub(".", "", annotations_full$Country, fixed = TRUE)
annotations_full$Country <- gsub(" ", "", annotations_full$Country, fixed = TRUE)
prob_nodes_full <- left_join(annotations_full, prob_long_full, by = c("node", "Country"))

#Let's remove all anotations that have less than 0.5 probability 
prob_nodes_0.5_full <- prob_nodes_full %>% filter(prob_nodes_full$value > 0.5)

#join with edge_table
edge_table_full <- left_join(edge_table_full, prob_nodes_0.5_full, by = 'node')
edge_table_full <- arrange(edge_table_full,child)


edge_table_full$tips<- sapply(edge_table_full$child, function(x){
  phangorn::Descendants(tree_full, x, type = 'tips')
})

node_full<-edge_table_full
for(i in 1:nrow(edge_table_full)) {
  node_full$key_lists[i] <- list(edge_table_full$node[unlist(edge_table_full$tips[i])])
}

names(node_full)[1] <- 'node_key'

predicted_full <- prob_nodes_0.5_full %>% filter(substr(node, 1, 1) == "n") %>% tally()
##Numbers to report:
paste("Out of", length(tree_full$node.label) , 'internal nodes in the full-canine tree, PastML is able to predict with 50% probability the ancestral character of',
        predicted_full$n,
        'nodes')

###################################
### Now repeat but with subsamples
###################################

## Sub1 

edge_table_sub1 <- data.frame(
  "child" = tree_sub1$edge[,2],
  "node" = sapply(tree_sub1$edge[,2],
                  select.tip.or.node,
                  tree = tree_sub1))

#Let's see all the node character predictions and the probabilities 
prob_long_sub1 <- melt(probabilities_sub1, id="node")
prob_long_sub1$value <- as.numeric(prob_long_sub1$value)
names(prob_long_sub1)[2] <- 'Country'
prob_long_sub1$Country <- gsub(".", "", prob_long_sub1$Country, fixed = TRUE)
prob_long_sub1$Country <- gsub(" ", "", prob_long_sub1$Country, fixed = TRUE)

#Let's join this with the annotations that have been attached for PastML
annotations_sub1$Country <- gsub(".", "", annotations_sub1$Country, fixed = TRUE)
annotations_sub1$Country <- gsub(" ", "", annotations_sub1$Country, fixed = TRUE)
prob_nodes_sub1 <- left_join(annotations_sub1, prob_long_sub1, by = c("node", "Country"))

#Let's remove all anotations that have less than 0.5 probability 
prob_nodes_0.5_sub1 <- prob_nodes_sub1 %>% filter(prob_nodes_sub1$value > 0.5)

#join with edge_table
edge_table_sub1 <- left_join(edge_table_sub1, prob_nodes_0.5_sub1, by = 'node')
edge_table_sub1 <- arrange(edge_table_sub1,child)


edge_table_sub1$tips<- sapply(edge_table_sub1$child, function(x){
  phangorn::Descendants(tree_sub1, x, type = 'tips')
})

node_sub1<-edge_table_sub1
for(i in 1:nrow(edge_table_sub1)) {
  #$key_lists[i] <- list(edge_table_sub1$node[grep(unlist(edge_table_sub1$tips[i]), edge_table_sub1$child)])
  node_sub1$key_lists[i] <- list(edge_table_sub1$node[unlist(edge_table_sub1$tips[i])])
}

predicted_sub1 <- prob_nodes_0.5_sub1 %>% filter(substr(node, 1, 1) == "n") %>% tally()
##Numbers to report:
paste("Out of", length(tree_sub1$node.label) , 'internal nodes in the sub1 tree, PastML is able to predict with 50% probability the ancestral character of',
      predicted_sub1$n,
      'nodes')


## Sub2 

edge_table_sub2 <- data.frame(
  "child" = tree_sub2$edge[,2],
  "node" = sapply(tree_sub2$edge[,2],
                  select.tip.or.node,
                  tree = tree_sub2))

#Let's see all the node character predictions and the probabilities 
prob_long_sub2 <- melt(probabilities_sub2, id="node")
prob_long_sub2$value <- as.numeric(prob_long_sub2$value)
names(prob_long_sub2)[2] <- 'Country'
prob_long_sub2$Country <- gsub(".", "", prob_long_sub2$Country, fixed = TRUE)
prob_long_sub2$Country <- gsub(" ", "", prob_long_sub2$Country, fixed = TRUE)

#Let's join this with the annotations that have been attached for PastML
annotations_sub2$Country <- gsub(".", "", annotations_sub2$Country, fixed = TRUE)
annotations_sub2$Country <- gsub(" ", "", annotations_sub2$Country, fixed = TRUE)
prob_nodes_sub2 <- left_join(annotations_sub2, prob_long_sub2, by = c("node", "Country"))

#Let's remove all anotations that have less than 0.5 probability 
prob_nodes_0.5_sub2 <- prob_nodes_sub2 %>% filter(prob_nodes_sub2$value > 0.5)

#join with edge_table
edge_table_sub2 <- left_join(edge_table_sub2, prob_nodes_0.5_sub2, by = 'node')
edge_table_sub2 <- arrange(edge_table_sub2,child)


edge_table_sub2$tips<- sapply(edge_table_sub2$child, function(x){
  phangorn::Descendants(tree_sub2, x, type = 'tips')
})

node_sub2<-edge_table_sub2
for(i in 1:nrow(edge_table_sub2)) {
  #$key_lists[i] <- list(edge_table_sub2$node[grep(unlist(edge_table_sub2$tips[i]), edge_table_sub2$child)])
  node_sub2$key_lists[i] <- list(edge_table_sub2$node[unlist(edge_table_sub2$tips[i])])
}

predicted_sub2 <- prob_nodes_0.5_sub2 %>% filter(substr(node, 1, 1) == "n") %>% tally()
##Numbers to report:
paste("Out of", length(tree_sub2$node.label) , 'internal nodes in the sub2 tree, PastML is able to predict with 50% probability the ancestral character of',
      predicted_sub2$n,
      'nodes')

## sub3 

edge_table_sub3 <- data.frame(
  "child" = tree_sub3$edge[,2],
  "node" = sapply(tree_sub3$edge[,2],
                  select.tip.or.node,
                  tree = tree_sub3))

#Let's see all the node character predictions and the probabilities 
prob_long_sub3 <- melt(probabilities_sub3, id="node")
prob_long_sub3$value <- as.numeric(prob_long_sub3$value)
names(prob_long_sub3)[2] <- 'Country'
prob_long_sub3$Country <- gsub(".", "", prob_long_sub3$Country, fixed = TRUE)
prob_long_sub3$Country <- gsub(" ", "", prob_long_sub3$Country, fixed = TRUE)

#Let's join this with the annotations that have been attached for PastML
annotations_sub3$Country <- gsub(".", "", annotations_sub3$Country, fixed = TRUE)
annotations_sub3$Country <- gsub(" ", "", annotations_sub3$Country, fixed = TRUE)
prob_nodes_sub3 <- left_join(annotations_sub3, prob_long_sub3, by = c("node", "Country"))

#Let's remove all anotations that have less than 0.5 probability 
prob_nodes_0.5_sub3 <- prob_nodes_sub3 %>% filter(prob_nodes_sub3$value > 0.5)

#join with edge_table
edge_table_sub3 <- left_join(edge_table_sub3, prob_nodes_0.5_sub3, by = 'node')
edge_table_sub3 <- arrange(edge_table_sub3,child)


edge_table_sub3$tips<- sapply(edge_table_sub3$child, function(x){
  phangorn::Descendants(tree_sub3, x, type = 'tips')
})

node_sub3<-edge_table_sub3
for(i in 1:nrow(edge_table_sub3)) {
  node_sub3$key_lists[i] <- list(edge_table_sub3$node[unlist(edge_table_sub3$tips[i])])
}

predicted_sub3 <- prob_nodes_0.5_sub3 %>% filter(substr(node, 1, 1) == "n") %>% tally()
##Numbers to report:
paste("Out of", length(tree_sub3$node.label) , 'internal nodes in the sub3 tree, PastML is able to predict with 50% probability the ancestral character of',
      predicted_sub3$n,
      'nodes')

## sub4 

edge_table_sub4 <- data.frame(
  "child" = tree_sub4$edge[,2],
  "node" = sapply(tree_sub4$edge[,2],
                  select.tip.or.node,
                  tree = tree_sub4))

#Let's see all the node character predictions and the probabilities 
prob_long_sub4 <- melt(probabilities_sub4, id="node")
prob_long_sub4$value <- as.numeric(prob_long_sub4$value)
names(prob_long_sub4)[2] <- 'Country'
prob_long_sub4$Country <- gsub(".", "", prob_long_sub4$Country, fixed = TRUE)
prob_long_sub4$Country <- gsub(" ", "", prob_long_sub4$Country, fixed = TRUE)

#Let's join this with the annotations that have been attached for PastML
annotations_sub4$Country <- gsub(".", "", annotations_sub4$Country, fixed = TRUE)
annotations_sub4$Country <- gsub(" ", "", annotations_sub4$Country, fixed = TRUE)
prob_nodes_sub4 <- left_join(annotations_sub4, prob_long_sub4, by = c("node", "Country"))

#Let's remove all anotations that have less than 0.5 probability 
prob_nodes_0.5_sub4 <- prob_nodes_sub4 %>% filter(prob_nodes_sub4$value > 0.5)

#join with edge_table
edge_table_sub4 <- left_join(edge_table_sub4, prob_nodes_0.5_sub4, by = 'node')
edge_table_sub4 <- arrange(edge_table_sub4,child)


edge_table_sub4$tips<- sapply(edge_table_sub4$child, function(x){
  phangorn::Descendants(tree_sub4, x, type = 'tips')
})

node_sub4<-edge_table_sub4
for(i in 1:nrow(edge_table_sub4)) {
  node_sub4$key_lists[i] <- list(edge_table_sub4$node[unlist(edge_table_sub4$tips[i])])
}

predicted_sub4 <- prob_nodes_0.5_sub4 %>% filter(substr(node, 1, 1) == "n") %>% tally()
##Numbers to report:
paste("Out of", length(tree_sub4$node.label) , 'internal nodes in the sub4 tree, PastML is able to predict with 50% probability the ancestral character of',
      predicted_sub4$n,
      'nodes')

## Sub5 

edge_table_sub5 <- data.frame(
  "child" = tree_sub5$edge[,2],
  "node" = sapply(tree_sub5$edge[,2],
                  select.tip.or.node,
                  tree = tree_sub5))

#Let's see all the node character predictions and the probabilities 
prob_long_sub5 <- melt(probabilities_sub5, id="node")
prob_long_sub5$value <- as.numeric(prob_long_sub5$value)
names(prob_long_sub5)[2] <- 'Country'
prob_long_sub5$Country <- gsub(".", "", prob_long_sub5$Country, fixed = TRUE)
prob_long_sub5$Country <- gsub(" ", "", prob_long_sub5$Country, fixed = TRUE)

#Let's join this with the annotations that have been attached for PastML
annotations_sub5$Country <- gsub(".", "", annotations_sub5$Country, fixed = TRUE)
annotations_sub5$Country <- gsub(" ", "", annotations_sub5$Country, fixed = TRUE)
prob_nodes_sub5 <- left_join(annotations_sub5, prob_long_sub5, by = c("node", "Country"))

#Let's remove all anotations that have less than 0.5 probability 
prob_nodes_0.5_sub5 <- prob_nodes_sub5 %>% filter(prob_nodes_sub5$value > 0.5)

#join with edge_table
edge_table_sub5 <- left_join(edge_table_sub5, prob_nodes_0.5_sub5, by = 'node')
edge_table_sub5 <- arrange(edge_table_sub5,child)


edge_table_sub5$tips<- sapply(edge_table_sub5$child, function(x){
  phangorn::Descendants(tree_sub5, x, type = 'tips')
})

node_sub5<-edge_table_sub5
for(i in 1:nrow(edge_table_sub5)) {
  node_sub5$key_lists[i] <- list(edge_table_sub5$node[unlist(edge_table_sub5$tips[i])])
}

predicted_sub5 <- prob_nodes_0.5_sub5 %>% filter(substr(node, 1, 1) == "n") %>% tally()
##Numbers to report:
paste("Out of", length(tree_sub5$node.label) , 'internal nodes in the sub5 tree, PastML is able to predict with 50% probability the ancestral character of',
      predicted_sub5$n,
      'nodes')

#########################################################################
# Now we want to compare across the trees and attach node_key identifier
#########################################################################

## Sub1 

setDT(node_full)[, l:=sapply(key_lists,length)]
f <- function(k) node_full[sapply(node_full$key_lists,\(i) all(k %chin% i))][l==min(l),node_key]
setDT(node_sub1)[, c:=sapply(key_lists,f)]
#recoding NA
node_sub1$Country <- ifelse(is.na(node_sub1$Country), 'unresolved', node_sub1$Country)
#identifying most frequent Country per node
node_sub1 <- node_sub1 %>% 
  group_by(c) %>% 
  summarise(Country = {t <- table(Country); names(t)[which.max(t)] })


## Sub2

setDT(node_full)[, l:=sapply(key_lists,length)]
f <- function(k) node_full[sapply(node_full$key_lists,\(i) all(k %chin% i))][l==min(l),node_key]
setDT(node_sub2)[, c:=sapply(key_lists,f)]
#recoding NA
node_sub2$Country <- ifelse(is.na(node_sub2$Country), 'unresolved', node_sub2$Country)
#identifying most frequent Country per node
node_sub2 <- node_sub2 %>% 
  group_by(c) %>% 
  summarise(Country = {t <- table(Country); names(t)[which.max(t)] })

## Sub3

setDT(node_full)[, l:=sapply(key_lists,length)]
f <- function(k) node_full[sapply(node_full$key_lists,\(i) all(k %chin% i))][l==min(l),node_key]
setDT(node_sub3)[, c:=sapply(key_lists,f)]
#recoding NA
node_sub3$Country <- ifelse(is.na(node_sub3$Country), 'unresolved', node_sub3$Country)
#identifying most frequent Country per node
node_sub3 <- node_sub3 %>% 
  group_by(c) %>% 
  summarise(Country = {t <- table(Country); names(t)[which.max(t)] })


## Sub4

setDT(node_full)[, l:=sapply(key_lists,length)]
f <- function(k) node_full[sapply(node_full$key_lists,\(i) all(k %chin% i))][l==min(l),node_key]
setDT(node_sub4)[, c:=sapply(key_lists,f)]
#recoding NA
node_sub4$Country <- ifelse(is.na(node_sub4$Country), 'unresolved', node_sub4$Country)
#identifying most frequent Country per node
node_sub4 <- node_sub4 %>% 
  group_by(c) %>% 
  summarise(Country = {t <- table(Country); names(t)[which.max(t)] })

## Sub5

setDT(node_full)[, l:=sapply(key_lists,length)]
f <- function(k) node_full[sapply(node_full$key_lists,\(i) all(k %chin% i))][l==min(l),node_key]
setDT(node_sub5)[, c:=sapply(key_lists,f)]
#recoding NA
node_sub5$Country <- ifelse(is.na(node_sub5$Country), 'unresolved', node_sub5$Country)
#identifying most frequent Country per node
node_sub5 <- node_sub5 %>% 
  group_by(c) %>% 
  summarise(Country = {t <- table(Country); names(t)[which.max(t)] })



# Now compare all subsamples and the full tree to get a consensus 

#remove duplicate keynodes from each subsample


names(node_sub1)[2] <- 'Country_Sub1'
names(node_sub2)[2] <- 'Country_Sub2'
names(node_sub3)[2] <- 'Country_Sub3'
names(node_sub4)[2] <- 'Country_Sub4'
names(node_sub5)[2] <- 'Country_Sub5'
names(node_full)[3] <- 'Country_Full'
names(node_full)[1] <- 'c'

node_sub1$c <- as.integer(node_sub1$c)
node_sub2$c <- as.integer(node_sub2$c)
node_sub3$c <- as.integer(node_sub3$c)
node_sub4$c <- as.integer(node_sub4$c)
node_sub5$c <- as.integer(node_sub5$c)


node_full_final <- node_full %>% left_join(node_sub1, by = 'c')
node_full_final <- node_full_final %>% left_join(node_sub2, by = 'c')
node_full_final <- node_full_final %>% left_join(node_sub3, by = 'c') 
node_full_final <- node_full_final %>% left_join(node_sub4, by = 'c')  
node_full_final <- node_full_final %>% left_join(node_sub5, by = 'c') 
  
node_full_final_final <- node_full_final %>% 
  select(c, l, node, Country_Full, Country_Sub1, Country_Sub2, Country_Sub3,
         Country_Sub4, Country_Sub5) %>% arrange(desc(l)) %>% filter(!is.na(Country_Full))


##########################################
## Some analysis now on this data
##########################################


node_comparison <- node_full_final_final
nrow(node_comparison)

node_comparison_NOLEAVE <- node_comparison %>% filter(l != 1)

nodes_inSub <- node_comparison %>% 
  filter(!is.na(Country_Sub1) | !is.na(Country_Sub2) | !is.na(Country_Sub3)
                | !is.na(Country_Sub4) | !is.na(Country_Sub5))

nodes_inSub_NOLEAVE <- nodes_inSub %>% filter(l != 1)

paste0(nrow(nodes_inSub_NOLEAVE)," out of ", nrow(node_comparison_NOLEAVE),
             " total nodes can be compared having found a similar node in subsampled trees")


nodes_sub_agg <- nodes_inSub %>% 
  select(c, Country_Full, Country_Sub1, Country_Sub2, Country_Sub3, Country_Sub4, Country_Sub5)

nodes_subOnly_agg <- nodes_sub_agg %>% 
  select(Country_Sub1, Country_Sub2, Country_Sub3, Country_Sub4, Country_Sub5)

agg <- apply(nodes_subOnly_agg,1,function(x) names(which.max(table(x))))
nodes_sub_agg$agg <- agg

nodes_inSub$agg <- nodes_sub_agg$agg 

nodes_inSub <- nodes_inSub %>% select(c,l,node, Country_Full, agg)

nodes_inSub$full_sub_same <- 
  ifelse(nodes_inSub$Country_Full == nodes_inSub$agg,
         TRUE, FALSE)
nodes_full_sub_same <- nodes_inSub %>% group_by(full_sub_same) %>% tally()
paste0("Out of ", nrow(nodes_inSub), " nodes ",
       nodes_full_sub_same$n[2], " nodes have the same ACR result between subsampled",
       " trees and the original full tree")

node_comparison %>% write.csv("/Volumes/NGS_Viroscreen/aholtz/euroME/GlobalRabies/total/alignment/better_alignment/complete_tree/complete_dating/canine/node_comparisonSUB.csv",
                              quote = FALSE, row.names = FALSE)
nodes_inSub %>% write.csv("/Volumes/NGS_Viroscreen/aholtz/euroME/GlobalRabies/total/alignment/better_alignment/complete_tree/complete_dating/canine/node_comparisonSUB_inSUB.csv",
                              quote = FALSE, row.names = FALSE)


###########
##BUILDING CONSENSUS TREE FROM TIPS THAT ARE UNION OF SUBSAMPLED TREES and Full Tree
## Using the same tree architecture from full tree, but using ACR from aggregated conseus subsamples

subtotal_tips <- append(tree_sub1$tip.label, tree_sub2$tip.label)
subtotal_tips <- subtotal_tips %>% 
  append(tree_sub3$tip.label) %>% 
  append(tree_sub4$tip.label) %>% 
  append(tree_sub5$tip.label)


subtotal_tips <- unique(subtotal_tips)

pruned_full_tree <- tree_full %>% keep.tip(subtotal_tips)

write.tree(pruned_full_tree, file = "/Volumes/NGS_Viroscreen/aholtz/euroME/GlobalRabies/total/alignment/better_alignment/complete_tree/complete_dating/canine/Fullpruned_subsampledTips.nwk")


#Prepare nodes_inSub for melting for PastML 
##
##n0	Brazil
##n0	Cambodia
##n0	Canada
##n0	Central African Republic

pastML_states <- nodes_inSub %>% select(node, Country_Full, agg)

##Sub1 Serbia-Poland fixes:

de<-data.frame(c("n99", "n9", "n88", "n4343", "n7", "n10", "n14611461", "n77", "n7", "n66", "n6", "n8", "n312312", "n311"), 
               c("Poland","Poland","Poland","Poland","Serbia","Germany","China", "Poland", "Poland", "Iraq", "Iraq", "Poland", "Iran", "Iraq"),
               c("Germany","Serbia","Serbia","Germany","Poland", "Germany", "China", "Serbia", "Serbia", "Serbia", "Serbia", "Serbia", "Iraq", "Iraq"))

names(de)<-c("node","Country_Full", "agg")

pastML_states<- rbind(pastML_states, de, fill =TRUE)
pastML_states_orig <- pastML_states
#

pastML_states <- pastML_states %>%  melt(id="node")
pastML_states <- pastML_states %>% select(-variable) %>% distinct() %>% filter(value != 'unresolved')
pastML_states$value <- as.character(pastML_states$value)



#write_delim(pastML_states, "/Volumes/NGS_Viroscreen/aholtz/euroME/GlobalRabies/total/alignment/better_alignment/complete_tree/complete_dating/canine/inSUB_Full_Sub_states.tab",
#            quote = 'none', col_names = TRUE, delim = "\t")
pastML_states <- read_delim("/Volumes/NGS_Viroscreen/aholtz/euroME/GlobalRabies/total/alignment/better_alignment/complete_tree/complete_dating/canine/inSUB_Full_Sub_states.tab")

pruned_full_tree <- read.tree("/Volumes/NGS_Viroscreen/aholtz/euroME/GlobalRabies/total/alignment/better_alignment/complete_tree/complete_dating/canine/Fullpruned_subsampledTips.nwk")

tips_missing<- as.data.frame(pruned_full_tree$tip.label) %>% subset(!pruned_full_tree$tip.label %in% pastML_states$node)
names(tips_missing)[1] <- "Accession"
tips_missing <- left_join(tips_missing, meta, by = "Accession")
tips_missing <- select(tips_missing, Accession, Country)

nodes_missing<- as.data.frame(pruned_full_tree$node.label) %>% subset(!(pruned_full_tree$node.label %in% pastML_states$node))
names(nodes_missing)[1] <- 'Accession'
nodes_missing$Country <- 'unresolved'

total_missing <- bind_rows(nodes_missing, tips_missing)
names(total_missing)[1] <- 'node'
names(total_missing)[2] <- 'value'

de_orig <-de
de <- de %>%  melt(id="node")
de <- select(de, -variable) %>% distinct()

pastML_states <- bind_rows(total_missing, pastML_states,de)

pastML_states$value <- gsub(" ", "", pastML_states$value)

write_delim(pastML_states, "/Volumes/NGS_Viroscreen/aholtz/euroME/GlobalRabies/total/alignment/better_alignment/complete_tree/complete_dating/canine/inSUB_Full_Sub_states.tab",
            quote = 'none', col_names = TRUE, delim = "\t")


##########
## Using the same tree but showing two columns
pastML_states_2Col <- pastML_states_orig
write_delim(pastML_states_2Col, "/Volumes/NGS_Viroscreen/aholtz/euroME/GlobalRabies/total/alignment/better_alignment/complete_tree/complete_dating/canine/inSUB_Full_Sub_states_2Col.tab",
            quote = 'none', col_names = TRUE, delim = "\t")
pastML_states_2Col <- read_delim("/Volumes/NGS_Viroscreen/aholtz/euroME/GlobalRabies/total/alignment/better_alignment/complete_tree/complete_dating/canine/inSUB_Full_Sub_states_2Col.tab")

tips_missing_2Col <- tips_missing
names(tips_missing_2Col)[1] <- 'node'
names(tips_missing_2Col)[2] <- 'Country_Full'
tips_missing_2Col$agg <- tips_missing_2Col$Country_Full

nodes_missing_2Col <- nodes_missing
names(nodes_missing_2Col)[1] <- 'node'
names(nodes_missing_2Col)[2] <- 'Country_Full'
nodes_missing_2Col$agg <- nodes_missing_2Col$Country_Full

pastML_states_2Col <- bind_rows(tips_missing_2Col,pastML_states_2Col,nodes_missing_2Col, de_orig)
pastML_states_2Col$Country_Full <- gsub(" ", "", pastML_states_2Col$Country_Full)
pastML_states_2Col$agg <- gsub(" ", "", pastML_states_2Col$agg)

write_delim(pastML_states_2Col, "/Volumes/NGS_Viroscreen/aholtz/euroME/GlobalRabies/total/alignment/better_alignment/complete_tree/complete_dating/canine/inSUB_Full_Sub_states_2Col.tab",
            quote = 'none', col_names = TRUE, delim = "\t")


