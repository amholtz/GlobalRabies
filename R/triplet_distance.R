####################################################
##                           Triplet Distances
## ####   Calculate Triplet Distances between subsampled tree 
## ####   by IQTREE + Gene Partitioning 
## ####   and 
## ####   Full Tree from FastTree reconstruction without partitioning
## author: Andrew Holtz
## creation date: 2022/03/16
###################################################



library(optparse)
library(Quartet)
library(tidyverse)
library(ape)

`%!in%` <- Negate(`%in%`)

##Parsing

option_list = list(
  make_option(c("-t", "--subtree"), type="character", default=NULL, 
              help="subtree NWK file path", metavar="character"),
  make_option(c("-f", "--fulltree"), type="character", default=NULL, 
              help="full bifurcated tree pruned for subtree tips NEXUS file path", metavar="character"),
  make_option(c("-n", "--tips"), type="double", default=NULL, 
              help="number of tips in subtree", metavar="double"),
  make_option(c("-o", "--table"), type="character", default=NULL, 
              help="out table file name", metavar="character")
);


opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

fulltree <- read.tree(opt$fulltree)
subtree <- read.tree(opt$subtree)

#Look for misplaced tips between trees

sub_tips <- subtree$tip.label
full_tips <- fulltree$tip.label
diff <- setdiff(sub_tips,full_tips)

new_sub <- drop.tip(subtree, diff)
write.tree(new_sub, 'temp_subTree.nwk')

# FullTree has been adapted to only contain tips in the subsample and is converted
# to a bifurcating tree

tDist <- TripletDistance(opt$fulltree, 'temp_subTree.nwk')

pND <- tDist/choose(opt$tips, 3)

tips <- opt$tips
binomial <- choose(opt$tips, 3)


out <- data.frame(tips,tDist,binomial,pND)

write.csv(out, opt$table,
            row.names = FALSE, quote = FALSE, col.names = TRUE)



