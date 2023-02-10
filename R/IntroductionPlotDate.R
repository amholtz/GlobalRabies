####################################################
##                           Timeline of Emergence
## ####   Comparing Emergence of Clades by Troupin, Villa, Holtz
## ####   
## author: Andrew Holtz
## creation date: 2023-JAN-27
###################################################

library(ggplot2)
library(tidyverse)
library(lubridate)
library(wesanderson)

setwd("/Volumes/NGS_Viroscreen/aholtz/euroME/project/GlobalRabies/")

intro <- read.csv(
  "/Volumes/NGS_Viroscreen/aholtz/euroME/project/GlobalRabies/data/SpecificFULLTREE_annotations_table_intro_nonneigh.csv")

intro$parent_date <- as.Date(date_decimal(intro$parent_date))
intro$child_date <- as.Date(date_decimal(intro$child_date))
intro$coloring <- paste(intro$parent_country, "->", intro$child_country)

intro=intro[order(intro$parent_date),]
intro$label=factor(intro$label,levels=intro$label)

ggplot(intro, aes(x = parent_date, y = label)) + 
  geom_linerange(aes(xmin = parent_date, xmax = child_date, size = subtree_size, color = speed)) +
  scale_colour_gradient(low = "#13005A", high = "#03C988", na.value = NA,limits = c(0,2000000))+
  scale_x_date(
    limits = c(as.Date("1890-01-01"), as.Date("2020-01-01")),
    date_break = "10 years", 
    date_labels = "%Y",
    expand = c(0,0)) +
  theme_linedraw() +
  xlab('Date Range of Transmission') +
  ylab('Country of Origin & Destination') +
  labs(size='Size of Subtree') +
  labs(color = 'km/year')


