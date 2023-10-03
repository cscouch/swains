#script to run two-factor anova in survey package on univariate metrics of colony size

#load libraries
library(tidyverse)
library(survey)
library(car) #Anova w/ Type III SS function
library(emmeans) #post-hoc options for Anova 
library(lme4)
library(multcomp)
library(hrbrthemes)
library(viridis)
library(ggridges)

rm(list=ls())
####Read in binned data and mean.SD file---------------
setwd("C:/github/swains/colony_size")
MOSP<-read.csv("MOSP_binned_tail_bins.csv", stringsAsFactors=FALSE) %>% mutate_if(is.character,as.factor)
POCS<-read.csv("POCS_binned_tail_bins.csv", stringsAsFactors=FALSE) %>% mutate_if(is.character,as.factor)
POSP<-read.csv("POSP_binned_tail_bins.csv", stringsAsFactors=FALSE) %>% mutate_if(is.character,as.factor)

MOSP$OBS_YEAR <- as.factor(MOSP$OBS_YEAR)
POCS$OBS_YEAR <- as.factor(POCS$OBS_YEAR)
POSP$OBS_YEAR <- as.factor(POSP$OBS_YEAR)

dat<-read.csv ("Adult_mean_SD.csv", stringsAsFactors=FALSE, row.names = NULL) %>% mutate_if(is.character,as.factor)
dat$OBS_YEAR <- as.factor(dat$OBS_YEAR)


#### Set survey design -------------------------------------------------------
#read in sector-area file
setwd("C:/github/swains/data")
sectors<-read.csv("Sectors-Strata-Areas.csv", stringsAsFactors=FALSE)
swa_sa<-filter(sectors,ISLAND=="Swains")

NH <- swa_sa %>%
  group_by(SEC_NAME, DEPTH_BIN)%>%
  summarize(unique = NH)%>%
  group_by(DEPTH_BIN)%>%
  summarise(NH = sum(unique))

MOSP<-left_join(MOSP,NH, by = "DEPTH_BIN")
POCS<-left_join(POCS,NH, by = "DEPTH_BIN")
POSP<-left_join(POSP,NH, by = "DEPTH_BIN")
dat<-left_join(dat,NH, by = "DEPTH_BIN")

#Calculate survey weights (inverse proportion weighting)
#MOSP
w.df<-MOSP %>%
  group_by(OBS_YEAR,REEF_ZONE,DEPTH_BIN,NH) %>%
  summarise(n = length(unique(SITE)))

w.df$sw<-w.df$NH/w.df$n #calculate survey weights for each site
MOSP.sw<-left_join(MOSP,w.df) #merge weights with site-level data
MOSP.sw$DEPTH_BIN <- as.factor(MOSP.sw$DEPTH_BIN)
head(MOSP.sw)
nrow(MOSP.sw) == n_distinct(MOSP.sw$SITE) #Check for duplicate sites

MOSP.sw$Strat_conc<-as.factor(paste(MOSP.sw$OBS_YEAR, MOSP.sw$DEPTH_BIN,sep = "_"))
#gather size classes into single response variable
MOSP.sw <- MOSP.sw %>%
  dplyr::select(DEPTH_BIN, OBS_YEAR, QJuv.R, Q10.R, QMed.R, Q90.R, sw, Strat_conc) %>%
  gather ("SIZE_BIN", "PROP", -DEPTH_BIN, -OBS_YEAR, -sw, -Strat_conc) %>% 
  mutate_if(is.character,as.factor)
MOSP.des<-svydesign(id=~1, strata=~ Strat_conc, weights=~sw,data=MOSP.sw)


#POCS
w.df<-POCS %>%
  group_by(OBS_YEAR,REEF_ZONE,DEPTH_BIN,NH) %>%
  summarise(n = length(unique(SITE)))

w.df$sw<-w.df$NH/w.df$n #calculate survey weights for each site
POCS.sw<-left_join(POCS,w.df) #merge weights with site-level data
POCS.sw$DEPTH_BIN <- as.factor(POCS.sw$DEPTH_BIN)
head(POCS.sw)
nrow(POCS.sw) == n_distinct(POCS.sw$SITE) #Check for duplicate sites

POCS.sw$Strat_conc<-as.factor(paste(POCS.sw$OBS_YEAR, POCS.sw$DEPTH_BIN,sep = "_"))
POCS.sw <- POCS.sw %>%
  dplyr::select(DEPTH_BIN, OBS_YEAR, QJuv.R, Q10.R, QMed.R, Q90.R, sw, Strat_conc) %>%
  gather ("SIZE_BIN", "PROP", -DEPTH_BIN, -OBS_YEAR, -sw, -Strat_conc) %>%
  mutate_if(is.character,as.factor)
POCS.des<-svydesign(id=~1, strata=~ Strat_conc, weights=~sw,data=POCS.sw)


#POSP
w.df<-POSP %>%
  group_by(OBS_YEAR,REEF_ZONE,DEPTH_BIN,NH) %>%
  summarise(n = length(unique(SITE)))

w.df$sw<-w.df$NH/w.df$n #calculate survey weights for each site
POSP.sw<-left_join(POSP,w.df) #merge weights with site-level data
POSP.sw$DEPTH_BIN <- as.factor(POSP.sw$DEPTH_BIN)
head(POSP.sw)
nrow(POSP.sw) == n_distinct(POSP.sw$SITE) #Check for duplicate sites

POSP.sw<-filter(POSP.sw, n!= "1") #remove strata that have less than 2 sites

POSP.sw$Strat_conc<-(paste(POSP.sw$OBS_YEAR, POSP.sw$DEPTH_BIN,sep = "_"))
POSP.sw <- POSP.sw %>%
  dplyr::select(DEPTH_BIN, OBS_YEAR, QJuv.R, Q10.R, QMed.R, Q90.R, sw, Strat_conc) %>%
  gather ("SIZE_BIN", "PROP", -DEPTH_BIN, -OBS_YEAR, -sw, -Strat_conc) %>%
  mutate_if(is.character,as.factor)
POSP.des<-svydesign(id=~1, strata=~ Strat_conc, weights=~sw,data=POSP.sw)




#### MEAN SIZE ANALYSIS: two-factor (YEAR *DEPTH) using site-level data for each depth strata-------------------------------

#MOSP-----
#MOSP Shallow

#test parametric assumptions
with(subset(MOSP.sw, DEPTH_BIN == "Shallow"), tapply((PROP), OBS_YEAR:SIZE_BIN, shapiro.test))
bartlett.test(PROP ~ Strat_conc, subset(MOSP.sw, DEPTH_BIN == "Shallow")) 

#run model
modR.shallow<-svyglm(PROP ~SIZE_BIN*OBS_YEAR, design=subset(MOSP.des, DEPTH_BIN == "Shallow"))
MOSP.shal <- Anova(modR.shallow, type = 3, test.statistic = "F") #significant inxn
MOSP.shal

#subset again by size bin
modR.shallow<-svyglm(PROP ~OBS_YEAR, design=subset(MOSP.des, DEPTH_BIN == "Shallow" | SIZE_BIN == "QJuv.R"))
Anova(modR.shallow, type = 3, test.statistic = "F") #NS; p =0.99 

modR.shallow<-svyglm(PROP ~OBS_YEAR, design=subset(MOSP.des, DEPTH_BIN == "Shallow" | SIZE_BIN == "Q10.R"))
Anova(modR.shallow, type = 3, test.statistic = "F") #NS; p =0.99 

modR.shallow<-svyglm(PROP ~OBS_YEAR, design=subset(MOSP.des, DEPTH_BIN == "Shallow" | SIZE_BIN == "QMed.R"))
Anova(modR.shallow, type = 3, test.statistic = "F") #NS; #NS; p =0.99 

modR.shallow<-svyglm(PROP ~OBS_YEAR, design=subset(MOSP.des, DEPTH_BIN == "Shallow" | SIZE_BIN == "Q90.R"))
Anova(modR.shallow, type = 3, test.statistic = "F") #NS; p =0.99 


#MOSP Mid

#test parametric assumptions
with(subset(MOSP.sw, DEPTH_BIN == "Mid"), tapply((PROP), OBS_YEAR:SIZE_BIN, shapiro.test))
bartlett.test(PROP ~ Strat_conc, subset(MOSP.sw, DEPTH_BIN == "Mid")) 

#run model
modR.Mid<-svyglm(PROP ~SIZE_BIN*OBS_YEAR, design=subset(MOSP.des, DEPTH_BIN == "Mid"))
MOSP.mid <- Anova(modR.Mid, type = 3, test.statistic = "F")
MOSP.mid #no significant interaction

#MOSP Deep

#test parametric assumptions
with(subset(MOSP.sw, DEPTH_BIN == "Deep"), tapply((PROP), OBS_YEAR:SIZE_BIN, shapiro.test))
bartlett.test(PROP ~ Strat_conc, subset(MOSP.sw, DEPTH_BIN == "Deep")) 

#run model
modR.Deep<-svyglm(PROP ~SIZE_BIN*OBS_YEAR, design=subset(MOSP.des, DEPTH_BIN == "Deep"))
MOSP.Deep <- Anova(modR.Deep, type = 3, test.statistic = "F") 
MOSP.Deep #no significant interaction


####POSP----
#POSP Shallow

#test parametric assumptions
with(subset(POSP.sw, DEPTH_BIN == "Shallow"), tapply((PROP), OBS_YEAR:SIZE_BIN, shapiro.test))
bartlett.test(PROP ~ Strat_conc, subset(POSP.sw, DEPTH_BIN == "Shallow")) 

#run model
modR.shallow<-svyglm(PROP ~SIZE_BIN*OBS_YEAR, design=subset(POSP.des, DEPTH_BIN == "Shallow"))
POSP.shal <- Anova(modR.shallow, type = 3, test.statistic = "F") #significant inxn
POSP.shal #no significant interaction

#POSP Mid

#test parametric assumptions
with(subset(POSP.sw, DEPTH_BIN == "Mid"), tapply((PROP), OBS_YEAR:SIZE_BIN, shapiro.test))
bartlett.test(PROP ~ Strat_conc, subset(POSP.sw, DEPTH_BIN == "Mid")) 

#run model
modR.Mid<-svyglm(PROP ~SIZE_BIN*OBS_YEAR, design=subset(POSP.des, DEPTH_BIN == "Mid"))
POSP.mid <- Anova(modR.Mid, type = 3, test.statistic = "F")
POSP.mid #no significant interaction

#POSP Deep

#test parametric assumptions
with(subset(POSP.sw, DEPTH_BIN == "Deep"), tapply((PROP), OBS_YEAR:SIZE_BIN, shapiro.test))
bartlett.test(PROP ~ Strat_conc, subset(POSP.sw, DEPTH_BIN == "Deep")) 

#run model
modR.Deep<-svyglm(PROP ~SIZE_BIN*OBS_YEAR, design=subset(POSP.des, DEPTH_BIN == "Deep"))
POSP.Deep <- Anova(modR.Deep, type = 3, test.statistic = "F") 
POSP.Deep #no significant interaction




####POCS----
#POCS Shallow

#test parametric assumptions
with(subset(POCS.sw, DEPTH_BIN == "Shallow"), tapply((PROP), OBS_YEAR:SIZE_BIN, shapiro.test))
bartlett.test(PROP ~ Strat_conc, subset(POCS.sw, DEPTH_BIN == "Shallow")) 

#run model
modR.shallow<-svyglm(PROP ~SIZE_BIN*OBS_YEAR, design=subset(POCS.des, DEPTH_BIN == "Shallow"))
POCS.shal <- Anova(modR.shallow, type = 3, test.statistic = "F") #significant inxn
POCS.shal


#subset again by size bin
modR.shallow<-svyglm(PROP ~OBS_YEAR, design=subset(POCS.des, DEPTH_BIN == "Shallow" | SIZE_BIN == "QJuv.R"))
Anova(modR.shallow, type = 3, test.statistic = "F") #NS; p =0.99 

modR.shallow<-svyglm(PROP ~OBS_YEAR, design=subset(POCS.des, DEPTH_BIN == "Shallow" | SIZE_BIN == "Q10.R"))
Anova(modR.shallow, type = 3, test.statistic = "F") #NS; p =0.99 

modR.shallow<-svyglm(PROP ~OBS_YEAR, design=subset(POCS.des, DEPTH_BIN == "Shallow" | SIZE_BIN == "QMed.R"))
Anova(modR.shallow, type = 3, test.statistic = "F") #NS; #NS; p =0.99 

modR.shallow<-svyglm(PROP ~OBS_YEAR, design=subset(POCS.des, DEPTH_BIN == "Shallow" | SIZE_BIN == "Q90.R"))
Anova(modR.shallow, type = 3, test.statistic = "F") #NS; p =0.99 



#POCS Mid

#test parametric assumptions
with(subset(POCS.sw, DEPTH_BIN == "Mid"), tapply((PROP), OBS_YEAR:SIZE_BIN, shapiro.test))
bartlett.test(PROP ~ Strat_conc, subset(POCS.sw, DEPTH_BIN == "Mid")) 

#run model
modR.Mid<-svyglm(PROP ~SIZE_BIN*OBS_YEAR, design=subset(POCS.des, DEPTH_BIN == "Mid"))
POCS.mid <- Anova(modR.Mid, type = 3, test.statistic = "F")
POCS.mid #no significant interaction

#POCS Deep

#test parametric assumptions
with(subset(POCS.sw, DEPTH_BIN == "Deep"), tapply((PROP), OBS_YEAR:SIZE_BIN, shapiro.test))
bartlett.test(PROP ~ Strat_conc, subset(POCS.sw, DEPTH_BIN == "Deep")) 

#run model
modR.Deep<-svyglm(PROP ~SIZE_BIN*OBS_YEAR, design=subset(POCS.des, DEPTH_BIN == "Deep"))
POCS.Deep <- Anova(modR.Deep, type = 3, test.statistic = "F") 
POCS.Deep #N.A. Not enough data to run fully crossed model

#subset again by size bin
modR.Deep<-svyglm(PROP ~OBS_YEAR, design=subset(POCS.des, DEPTH_BIN == "Deep" | SIZE_BIN == "QJuv.R"))
Anova(modR.Deep, type = 3, test.statistic = "F") #NS; p =0.91

modR.Deep<-svyglm(PROP ~OBS_YEAR, design=subset(POCS.des, DEPTH_BIN == "Deep" | SIZE_BIN == "Q10.R"))
Anova(modR.Deep, type = 3, test.statistic = "F") #NS; p =0.052

modR.Deep<-svyglm(PROP ~OBS_YEAR, design=subset(POCS.des, DEPTH_BIN == "Deep" | SIZE_BIN == "QMed.R"))
Anova(modR.Deep, type = 3, test.statistic = "F") #NS; #NS; p =0.378 

modR.Deep<-svyglm(PROP ~OBS_YEAR, design=subset(POCS.des, DEPTH_BIN == "Deep" | SIZE_BIN == "Q90.R"))
Anova(modR.Deep, type = 3, test.statistic = "F") #NS; p =0.3145 









