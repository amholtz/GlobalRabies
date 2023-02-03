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

comp <- read.csv(
  "/Volumes/NGS_Viroscreen/aholtz/euroME/project/GlobalRabies/data/CITree_Comparison.csv")

comp$End<-lubridate::ymd(comp$End, truncated = 2L)
comp$Start<-lubridate::ymd(comp$Start, truncated = 2L)
comp$Point<-lubridate::ymd(comp$Point, truncated = 2L)

comp=comp[order(comp$present),]
comp$Node=factor(comp$Node,levels=comp$Node)

plot <- ggplot(comp, aes(x = Point, y = Node)) + 
  geom_linerange(aes(xmin = Start, xmax = End, color = paper), size = 2) +
  scale_colour_manual(values =wes_palette("Cavalcanti1")) + 
  geom_point(size = 2, fill = 'black', color = 'white') + 
  scale_x_date(
    limits = c(as.Date("1250-01-01"), as.Date("1880-01-01")),
    date_break = "100 years", 
    date_labels = "%Y",
    expand = c(0,0)) +
  theme_linedraw() +
  xlab('Year') +
  ylab('Study & Node') +
  labs(color='Study')

plot





