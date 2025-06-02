rm(list = ls())

library(gdata)             # needed for drop_levels()
library(reshape)           # reshape library inclues the cast() function used below
library(splitstackshape)
library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(survey)
library(multcomp)
library(emmeans)

dir = Sys.info()[7]
setwd(paste0("C:/Users/", dir, "/Documents/GitHub/swains/"))


#LOAD site-level data & Prep
swains_t2_pooled <- read.csv("Data/BenthicCover_2010-2023_Tier2b.csv")
swa<- swains_t2_pooled %>%
  dplyr::select(c(ISLAND,SITE,OBS_YEAR, new_MAX_DEPTH_M,DEPTH_BIN,POCS,POSP,MOSP)) %>%
  pivot_longer(cols= POCS:MOSP, names_to = "New",values_to = "Percent") %>%
  filter(ISLAND =="Swains" & new_MAX_DEPTH_M >=3) #subset sites deeper than 3m
head(swa)
nrow(swa)

#read in sector-area file
sectors<-read.csv("Data/Sectors-Strata-Areas.csv", stringsAsFactors=FALSE)
swa_sa<-filter(sectors,ISLAND=="Swains") #subset to swains only
swa_sa$DEPTH_BIN<-as.factor(swa_sa$DEPTH_BIN)

NH <- swa_sa %>%
  group_by(SEC_NAME, DEPTH_BIN)%>%
  dplyr::summarize(unique = NH)%>%
  group_by(DEPTH_BIN)%>%
  dplyr::summarise(NH = sum(unique)) #concatenate swains open and sanctuary sectors bc we don't care about differences between these areas

swa<-left_join(swa,NH) #merge percent cover data with new NH values pooled across the 2 swains sectors


#Calculate survey weights (inverse proportion weighting)
w.df<-swa %>%
  group_by(OBS_YEAR,DEPTH_BIN,NH) %>%
  dplyr::summarise(n = length(unique(SITE)))

w.df$sw<-w.df$NH/w.df$n #calculate survey weights for each site

site.sw<-left_join(swa,w.df) #merge weights with site-level data
head(site.sw)

#site.sw$Strat_conc<-paste(site.sw$OBS_YEAR, site.sw$DEPTH_BIN,sep = "_")
site.sw$OBS_YEAR<-as.character(site.sw$OBS_YEAR)
site.sw$OBS_YEAR<-as.factor(site.sw$OBS_YEAR)
site.sw$DEPTH_BIN<-as.factor(site.sw$DEPTH_BIN)

#create strata levels for survey design
site.sw.pool <- site.sw
site.sw.pool$Strat_conc<-paste(site.sw.pool$OBS_YEAR, site.sw.pool$DEPTH_BIN, site.sw.pool$New, sep = "_")#create strata design including taxa so we can survey lonely.psu correctly
options(survey.lonely.psu = "adjust") #keep in strata with only 1 site

des<-svydesign(id=~1, strata=~ Strat_conc, weights=~sw,data= filter(site.sw.pool, New!= "Other")) 

##########ALL TAXA BENTHIC COVER PLOT#####
#create strata level estimates of percent cover from the survey package
strata.means <- svyby(~Percent, ~Strat_conc, des, svymean) #calculate mean and standard error based on weighting
strata.means <- strata.means %>% tidyr::separate(Strat_conc, sep = "_", into = c("OBS_YEAR", "DEPTH_BIN", "TAXA")) #unconcatenate strata variables

library(forcats)
strata.means <- strata.means %>% 
  filter(OBS_YEAR %in% c("2015","2018", "2023")) %>%
  mutate(DEPTH_BIN=factor(DEPTH_BIN)) %>% 
  mutate(DEPTH_BIN=fct_relevel(DEPTH_BIN,c("Shallow", "Mid", "Deep"))) %>%
  arrange(DEPTH_BIN)


strata.means <- strata.means %>% 
  mutate(OBS_YEAR=factor(OBS_YEAR)) %>% 
  mutate(OBS_YEAR=fct_relevel(OBS_YEAR,c("2015","2018", "2023"))) %>%
  arrange(OBS_YEAR)


strata.means <- strata.means %>% 
  mutate(TAXA=factor(TAXA)) %>% 
  mutate(TAXA=fct_relevel(TAXA,c("POCS", "MOSP", "POSP"))) %>%
  arrange(TAXA)

strata.means$Percent<-round(strata.means$Percent,digits = 2)
strata.means$se<-round(strata.means$se,digits = 1)
strata.means$Mean_SE <- paste(strata.means$Percent,strata.means$se,sep= " ± ")

strata.means
write.csv(strata.means,"Data/Genus_levelStrata_level_benthic_cover.csv")
