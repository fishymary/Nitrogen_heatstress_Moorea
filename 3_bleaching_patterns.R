# -----------------------------------------------------------------------
# Nitrogen pollution interacts with thermal stress to increase coral bleaching across the seascape
# Donovan et al. 
# 3 - Bleaching patterns summary
# -----------------------------------------------------------------------

# This script summarizes bleaching patterns and creates rasters for mapping in Figure 3.

# Initialization ----------------------------------------------------------
library(dplyr)
library(rjags)
library(boot)

# data --------------------------------------------------------------------
# bleaching
moorea <- read.csv('data/moorea_withavgnuts.csv')

# summarize ---------------------------------------------------------------
poc_prev <- moorea %>% filter(Taxa=="Pocillopora" & !is.na(Percent_bleached)) %>% group_by(Point,Latitude,Longitude) %>% summarise('meanPrev'=mean(y)) %>% ungroup()
hist(poc_prev$meanPrev)
summary(poc_prev$meanPrev)
moorea %>% filter(Taxa=="Pocillopora" & !is.na(Percent_bleached)) %>% group_by(Island_shore,Point,Latitude,Longitude) %>% summarise('meanPrev'=mean(y)) %>% ungroup() %>% group_by(Island_shore) %>% summarise('mean'=mean(meanPrev),'max'=max(meanPrev),'min'=min(meanPrev)) %>% ungroup()

acr_prev <- moorea %>% filter(Taxa=="Acropora" & !is.na(Percent_bleached)) %>% group_by(Point,Latitude,Longitude) %>% summarise('meanPrev'=mean(y)) %>% ungroup()
hist(acr_prev$meanPrev)

poc_sev <- moorea %>% filter(Taxa=="Pocillopora" & !is.na(Percent_bleached) & y==1) %>% group_by(Point,Latitude,Longitude) %>% summarise('meanSev'=mean(Percent_bleached/100)) %>% ungroup()
summary(poc_sev$meanSev)
moorea %>% filter(Taxa=="Pocillopora" & !is.na(Percent_bleached) & y==1) %>% group_by(Island_shore,Point,Latitude,Longitude) %>% summarise('meanSev'=mean(Percent_bleached/100)) %>% ungroup() %>% group_by(Island_shore) %>% summarise('median'=median(meanSev),'max'=max(meanSev),'min'=min(meanSev)) %>% ungroup()

acr_sev <- moorea %>% filter(Taxa=="Acropora" & !is.na(Percent_bleached) & y==1) %>% group_by(Point,Latitude,Longitude) %>% summarise('meanSev'=mean(Percent_bleached/100)) %>% ungroup()
hist(acr_sev$meanSev)
summary(acr_sev$meanSev)
moorea %>% filter(Taxa=="Acropora" & !is.na(Percent_bleached) & y==1) %>% group_by(Island_shore,Point,Latitude,Longitude) %>% summarise('meanSev'=mean(Percent_bleached/100)) %>% ungroup() %>% group_by(Island_shore) %>% summarise('median'=median(meanSev),'max'=max(meanSev),'min'=min(meanSev)) %>% ungroup()
