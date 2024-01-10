rm(list = ls())

library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)

source("C:/Users/jonathan.charendoff/Documents/GitHub/Benthic-Scripts/Functions/Benthic_Functions_newApp.R")
source("C:/Users/jonathan.charendoff/Documents/GitHub/fish-paste/lib/core_functions.R")
source("C:/Users/jonathan.charendoff/Documents/GitHub/fish-paste/lib/fish_team_functions.R")
source("C:/Users/jonathan.charendoff/Documents/GitHub/fish-paste/lib/Islandwide Mean&Variance Functions.R")

#FUNCTIONAL GROUP LOOKUP TABLE
new.grouping <- read.csv("C:/Users/Jonathan.Charendoff/Documents/GitHub/swains/Scripts/Benthic Cover/Benthic_Cover_Lookup_Table.csv")

#read in site level data for each tier
t2 <- read.csv("T:/Benthic/Data/REA Coral Demography & Cover/Summary Data/Site/BenthicCover_2010-2023_Tier2b_SITE.csv")

#filter site data to only be swains, StRS and filter out columns we don't want
swains_t2_w <- t2[t2$ISLAND == "Swains" & t2$Oceanography == 0,-(106:109)]

#wide to long format
swains_t2 <- pivot_longer(swains_t2_w, cols = c(24:ncol(swains_t2_w)), names_to = "TAXA", values_to = "PERCENT")

#remove taxa that are never present 
test <- swains_t2 %>% 
  group_by(TAXA) %>% 
  dplyr::summarise(zeros = sum(PERCENT))
zeros <- test$TAXA[test$zeros == 0]
swains_t2 <- swains_t2[-(which(swains_t2$TAXA %in% zeros)),]

#update montastrea to astrea
swains_t2$TAXA <- recode(swains_t2$TAXA, "MONS" = "ASTS")

#NEW FUNCTIONAL GROUP LABELING
swains_t2 <- left_join(swains_t2, new.grouping) #merge in new labels
swains_t2$New[is.na(swains_t2$New)] <- "Hard Coral"
swains_t2$New <- recode_factor(swains_t2$New, "POSP"="Total Coral", "MOSP" = "Total Coral", "POCS" = "Total Coral",  
                               "CCA" = "CCA", "Hard Coral" = "Total Coral",
                               .default = "Other", "MICR" = "UPMA","UPMA" = "UPMA", "TURF" = "TURF", "EMA" = "EMA") #relabel the taxa codes to be more general, but can change if you want to looks at focal corals


#pool data with new functional grouping
swains_t2_pooled <- swains_t2 %>% 
  group_by(SITE, OBS_YEAR, DATE_, ANALYSIS_YEAR, new_MIN_DEPTH_M, new_MAX_DEPTH_M, LATITUDE, LONGITUDE, DEPTH_BIN, New) %>% 
  dplyr::summarise(Percent = sum(PERCENT))

#export cleaned data
write.csv(swains_t2_pooled, "C:/Users/Jonathan.Charendoff/Documents/GitHub/swains/Data/Benthic_Cover_Ready.csv")
