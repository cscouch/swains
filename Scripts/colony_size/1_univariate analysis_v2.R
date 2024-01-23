#script to run two-factor anova in survey package on univariate metrics of colony size

#load libraries
library(tidyverse)
library(survey)
library(car) #Anova w/ Type III SS function
library(emmeans) #post-hoc options for Anova 
library(lme4)
library(multcomp)

rm(list=ls())
####Read in binned data and mean.SD file---------------

dir = Sys.info()[7]
setwd(paste0("C:/Users/", dir, "/Documents/GitHub/swains/"))
#setwd("C:/github/swains") #BH directory

#read in data
MOSP<-read.csv("Data/MOSP_binned_tail_bins.csv", stringsAsFactors=FALSE) %>% mutate_if(is.character,as.factor)
POCS<-read.csv("Data/POCS_binned_tail_bins.csv", stringsAsFactors=FALSE) %>% mutate_if(is.character,as.factor)
POSP<-read.csv("Data/POSP_binned_tail_bins.csv", stringsAsFactors=FALSE) %>% mutate_if(is.character,as.factor)
sectors<-read.csv("Data/Sectors-Strata-Areas.csv", stringsAsFactors=FALSE)

MOSP$OBS_YEAR <- as.factor(MOSP$OBS_YEAR)
POCS$OBS_YEAR <- as.factor(POCS$OBS_YEAR)
POSP$OBS_YEAR <- as.factor(POSP$OBS_YEAR)


#### Set survey design -------------------------------------------------------
swa_sa<-filter(sectors,ISLAND=="Swains")

NH <- swa_sa %>%
  group_by(SEC_NAME, DEPTH_BIN)%>%
  summarize(unique = NH)%>%
  group_by(DEPTH_BIN)%>%
  summarise(NH = sum(unique))

MOSP<-left_join(MOSP,NH, by = "DEPTH_BIN")
POCS<-left_join(POCS,NH, by = "DEPTH_BIN")
POSP<-left_join(POSP,NH, by = "DEPTH_BIN")

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



#### 2-factor svyglm analysis (SIZE BIN * OBS_YEAR) for each depth strata and for each taxa-------------------------------

#MOSP-----
#MOSP Shallow
#test parametric assumptions
with(filter(MOSP.sw, DEPTH_BIN == "Shallow"), tapply((PROP), OBS_YEAR:SIZE_BIN, shapiro.test))
bartlett.test(PROP ~ Strat_conc, filter(MOSP.sw, DEPTH_BIN == "Shallow")) 

#run model
MOSP.sw_shal<-filter(MOSP.sw, DEPTH_BIN == "Shallow")
MOSP.des_shal<-svydesign(id=~1, strata=~ Strat_conc, weights=~sw,data=MOSP.sw_shal)
modR.shallow<-svyglm(PROP ~SIZE_BIN*OBS_YEAR, design=MOSP.des_shal) #default is gaussian family
Anova(modR.shallow, type = 3, test.statistic = "F") #significant inxn

#subset again by size bin
MOSP.sw_shal_juv<-filter(MOSP.sw_shal, SIZE_BIN == "QJuv.R")
MOSP.des_shal_juv<-svydesign(id=~1, strata=~ Strat_conc, weights=~sw,data=MOSP.sw_shal_juv)
modR.shallow<-svyglm(PROP ~OBS_YEAR, design=MOSP.des_shal_juv)
Anova(modR.shallow, type = 3, test.statistic = "F") #NS; p =0.84 

MOSP.sw_shal_Q10<-filter(MOSP.sw_shal, SIZE_BIN == "Q10.R")
MOSP.des_shal_Q10<-svydesign(id=~1, strata=~ Strat_conc, weights=~sw,data=MOSP.sw_shal_Q10)
modR.shallow<-svyglm(PROP ~OBS_YEAR,  design=MOSP.des_shal_Q10)
Anova(modR.shallow, type = 3, test.statistic = "F") #p =0.001
emmeans(modR.shallow, pairwise ~ OBS_YEAR) #2015 < 2023

dat_weighted_MOSP<-svyby(~PROP,~OBS_YEAR,MOSP.des_shal_Q10,svymean)


MOSP.sw_shal_QMed<-filter(MOSP.sw_shal, SIZE_BIN == "QMed.R")
MOSP.des_shal_QMed<-svydesign(id=~1, strata=~ Strat_conc, weights=~sw,data=MOSP.sw_shal_QMed)
modR.shallow<-svyglm(PROP ~OBS_YEAR,  design=MOSP.des_shal_QMed)
Anova(modR.shallow, type = 3, test.statistic = "F") #NS;

MOSP.sw_shal_Q90<-filter(MOSP.sw_shal, SIZE_BIN == "Q90.R")
MOSP.des_shal_Q90<-svydesign(id=~1, strata=~ Strat_conc, weights=~sw,data=MOSP.sw_shal_Q90)
modR.shallow<-svyglm(PROP ~OBS_YEAR,  design=MOSP.des_shal_Q90)
Anova(modR.shallow, type = 3, test.statistic = "F") #NS;


#MOSP Mid
#test parametric assumptions
with(filter(MOSP.sw, DEPTH_BIN == "Mid"), tapply((PROP), OBS_YEAR:SIZE_BIN, shapiro.test))
bartlett.test(PROP ~ Strat_conc, filter(MOSP.sw, DEPTH_BIN == "Mid")) 

#run model
MOSP.sw_mid<-filter(MOSP.sw, DEPTH_BIN == "Mid")
MOSP.des_mid<-svydesign(id=~1, strata=~ Strat_conc, weights=~sw,data=MOSP.sw_mid)
modR.mid<-svyglm(PROP ~SIZE_BIN*OBS_YEAR, design=MOSP.des_mid) 
Anova(modR.mid, type = 3, test.statistic = "F") #NS inxn


#MOSP Deep
#test parametric assumptions
with(subset(MOSP.sw, DEPTH_BIN == "Deep"), tapply((PROP), OBS_YEAR:SIZE_BIN, shapiro.test))
bartlett.test(PROP ~ Strat_conc, subset(MOSP.sw, DEPTH_BIN == "Deep")) 

#run model
MOSP.sw_deep<-filter(MOSP.sw, DEPTH_BIN == "Deep")
MOSP.des_deep<-svydesign(id=~1, strata=~ Strat_conc, weights=~sw,data=MOSP.sw_deep)
modR.deep<-svyglm(PROP ~SIZE_BIN*OBS_YEAR, design=MOSP.des_deep) 
Anova(modR.deep, type = 3, test.statistic = "F") #NS inxn; Year is significant
emmeans(modR.deep, pairwise ~ OBS_YEAR) #



####POSP----
#POSP Shallow
#test parametric assumptions
with(subset(POSP.sw, DEPTH_BIN == "Shallow"), tapply((PROP), OBS_YEAR:SIZE_BIN, shapiro.test))
bartlett.test(PROP ~ Strat_conc, subset(POSP.sw, DEPTH_BIN == "Shallow")) 

#run model
POSP.sw_shal<-filter(POSP.sw, DEPTH_BIN == "Shallow")
POSP.des_shal<-svydesign(id=~1, strata=~ Strat_conc, weights=~sw,data=POSP.sw_shal)
modR.shallow<-svyglm(PROP ~SIZE_BIN*OBS_YEAR, design=POSP.des_shal)
Anova(modR.shallow, type = 3, test.statistic = "F") #NS inxn


#POSP Mid
#test parametric assumptions
with(subset(POSP.sw, DEPTH_BIN == "Mid"), tapply((PROP), OBS_YEAR:SIZE_BIN, shapiro.test))
bartlett.test(PROP ~ Strat_conc, subset(POSP.sw, DEPTH_BIN == "Mid")) 

#run model
POCS.sw_mid<-filter(POSP.sw, DEPTH_BIN == "Mid")
POSP.des_mid<-svydesign(id=~1, strata=~ Strat_conc, weights=~sw,data=POSP.sw_mid)
modR.mid<-svyglm(PROP ~SIZE_BIN*OBS_YEAR, design=POSP.des_mid)
Anova(modR.mid, type = 3, test.statistic = "F") #NS inxn


#POSP Deep
#test parametric assumptions
with(subset(POSP.sw, DEPTH_BIN == "Deep"), tapply((PROP), OBS_YEAR:SIZE_BIN, shapiro.test))
bartlett.test(PROP ~ Strat_conc, subset(POSP.sw, DEPTH_BIN == "Deep")) 

#run model
POSP.sw_deep<-filter(POSP.sw, DEPTH_BIN == "Deep")
POSP.des_deep<-svydesign(id=~1, strata=~ Strat_conc, weights=~sw,data=POSP.sw_deep)
modR.deep<-svyglm(PROP ~SIZE_BIN*OBS_YEAR, design=POSP.des_deep) 
Anova(modR.deep, type = 3, test.statistic = "F")  #NS inxn



####POCS----
#POCS Shallow

#test parametric assumptions
with(subset(POCS.sw, DEPTH_BIN == "Shallow"), tapply((PROP), OBS_YEAR:SIZE_BIN, shapiro.test))
bartlett.test(PROP ~ Strat_conc, subset(POCS.sw, DEPTH_BIN == "Shallow")) 

#run model
POCS.sw_shal<-filter(POCS.sw, DEPTH_BIN == "Shallow")
POCS.des_shal<-svydesign(id=~1, strata=~ Strat_conc, weights=~sw,data=POCS.sw_shal)
modR.shallow<-svyglm(PROP ~SIZE_BIN*OBS_YEAR, design=POCS.des_shal) 
Anova(modR.shallow, type = 3, test.statistic = "F") #significant inxn

#subset again by size bin
POCS.sw_shal_juv<-filter(POCS.sw_shal, SIZE_BIN == "QJuv.R")
POCS.des_shal_juv<-svydesign(id=~1, strata=~ Strat_conc, weights=~sw,data=POCS.sw_shal_juv)
modR.shallow<-svyglm(PROP ~OBS_YEAR, design=POCS.des_shal_juv)
Anova(modR.shallow, type = 3, test.statistic = "F") #NS 

POCS.sw_shal_Q10<-filter(POCS.sw_shal, SIZE_BIN == "Q10.R")
POCS.des_shal_Q10<-svydesign(id=~1, strata=~ Strat_conc, weights=~sw,data=POCS.sw_shal_Q10)
modR.shallow<-svyglm(PROP ~OBS_YEAR,  design=POCS.des_shal_Q10)
Anova(modR.shallow, type = 3, test.statistic = "F") #p =0.001
emmeans(modR.shallow, pairwise ~ OBS_YEAR) #2015 < 2023

dat_weighted<-svyby(~PROP,~OBS_YEAR,POCS.des_shal_Q10,svymean)

POCS.sw_shal_QMed<-filter(POCS.sw_shal, SIZE_BIN == "QMed.R")
POCS.des_shal_QMed<-svydesign(id=~1, strata=~ Strat_conc, weights=~sw,data=POCS.sw_shal_QMed)
modR.shallow<-svyglm(PROP ~OBS_YEAR,  design=POCS.des_shal_QMed)
Anova(modR.shallow, type = 3, test.statistic = "F") #NS

POCS.sw_shal_Q90<-filter(POCS.sw_shal, SIZE_BIN == "Q90.R")
POCS.des_shal_Q90<-svydesign(id=~1, strata=~ Strat_conc, weights=~sw,data=POCS.sw_shal_Q90)
modR.shallow<-svyglm(PROP ~OBS_YEAR,  design=POCS.des_shal_Q90)
Anova(modR.shallow, type = 3, test.statistic = "F") #p =0.001

#POCS Mid
#test parametric assumptions
with(subset(POCS.sw, DEPTH_BIN == "Mid"), tapply((PROP), OBS_YEAR:SIZE_BIN, shapiro.test))
bartlett.test(PROP ~ Strat_conc, subset(POCS.sw, DEPTH_BIN == "Mid")) 

#run model
POCS.sw_mid<-filter(POCS.sw, DEPTH_BIN == "Mid")
POCS.des_mid<-svydesign(id=~1, strata=~ Strat_conc, weights=~sw,data=POCS.sw_mid)
modR.mid<-svyglm(PROP ~SIZE_BIN*OBS_YEAR, design=POCS.des_mid) #default is gaussian famil
Anova(modR.mid, type = 3, test.statistic = "F") #NS


#POCS Deep
#test parametric assumptions
with(subset(POCS.sw, DEPTH_BIN == "Deep"), tapply((PROP), OBS_YEAR:SIZE_BIN, shapiro.test))
bartlett.test(PROP ~ Strat_conc, subset(POCS.sw, DEPTH_BIN == "Deep")) 

#run model
POCS.sw_deep<-filter(POCS.sw, DEPTH_BIN == "Deep")
POCS.des_deep<-svydesign(id=~1, strata=~ Strat_conc, weights=~sw,data=POCS.sw_deep)
modR.deep<-svyglm(PROP ~SIZE_BIN*OBS_YEAR, design=POCS.des_deep) #default is gaussian famil
Anova(modR.deep, type = 3, test.statistic = "F") #NS
#N.A. Not enough data to run fully crossed model

#subset again by size bin
POCS.sw_deep_juv<-filter(POCS.sw_deep, SIZE_BIN == "QJuv.R")
POCS.des_deep_juv<-svydesign(id=~1, strata=~ Strat_conc, weights=~sw,data=POCS.sw_deep_juv)
svyranktest(PROP ~ OBS_YEAR, design=POCS.des_deep_juv, test=("KruskalWallis"))
modR.deep<-svyglm(PROP ~OBS_YEAR, design=POCS.des_deep_juv)
Anova(modR.deep, type = 3, test.statistic = "F") #NS 
#both the parametric and non-parametric approaches are throwing singularity errors; lack o' data (aka only one site with density of POCS)

POCS.sw_deep_Q10<-filter(POCS.sw_deep, SIZE_BIN == "Q10.R")
POCS.des_deep_Q10<-svydesign(id=~1, strata=~ Strat_conc, weights=~sw,data=POCS.sw_deep_Q10)
modR.deep<-svyglm(PROP ~OBS_YEAR,  design=POCS.des_deep_Q10)
Anova(modR.deep, type = 3, test.statistic = "F") #NS

POCS.sw_deep_QMed<-filter(POCS.sw_deep, SIZE_BIN == "QMed.R")
POCS.des_deep_QMed<-svydesign(id=~1, strata=~ Strat_conc, weights=~sw,data=POCS.sw_deep_QMed)
modR.deep<-svyglm(PROP ~OBS_YEAR,  design=POCS.des_deep_QMed)
Anova(modR.deep, type = 3, test.statistic = "F") #NS

POCS.sw_deep_Q90<-filter(POCS.sw_deep, SIZE_BIN == "Q90.R")
POCS.des_deep_Q90<-svydesign(id=~1, strata=~ Strat_conc, weights=~sw,data=POCS.sw_deep_Q90)
svyranktest(PROP ~ OBS_YEAR, design=POCS.des_deep_Q90, test=("KruskalWallis"))
modR.deep<-svyglm(PROP ~OBS_YEAR,  design=POCS.des_deep_Q90)
Anova(modR.deep, type = 3, test.statistic = "F") #p =0.001
#both the parametric and non-parametric approaches are throwing singularity errors; lack o' data (aka only one site with density of POCS)



