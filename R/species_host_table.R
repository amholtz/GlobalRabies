####################################################
##        Table with Host Species and Cleaned Host 
## ####   
## author: Andrew Holtz
## creation date: 2023/02/20
###################################################

rm(list=ls())
library(dplyr)
library(data.table)
library(taxize)
library(readr)
library(optparse)

##Parsing

option_list = list(
  make_option(c("-d", "--meta"), type="character", default=NULL, 
              help="metadata file path", metavar="character"),
  make_option(c("-o", "--table"), type="character", default=NULL, 
              help="output table path", metavar="character")
);

  opt_parser = OptionParser(option_list=option_list)
  opt = parse_args(opt_parser)
meta <- read.delim(opt$meta)
print(meta$Host)

host_species <- unique(meta$Host)

x <- taxize::tax_name(host_species, get = 'family', db = 'ncbi')
names(family)[2] <- 'species'
family <- family[-1]
family$order <- taxize::tax_name(family$family, get = 'order', db = 'ncbi')

family$keep <- ifelse(family$species == 'Canis lupus familiaris', 'dog',
                      ifelse(family$species == 'Canis lupus', 'wolf',
                             ifelse(family$species == 'Vulpes vulpes', 'fox',
                                    ifelse(family$order$order == 'Chiroptera', 'bat',
                                           family$family))))

meta$host_simple2 <- ifelse(meta$Host %in% family$species, family$keep, "other")

host <- meta %>% select(Host, host_simple2) %>% distinct()

write.csv(host, opt$table, row.names = FALSE, quote = FALSE)

