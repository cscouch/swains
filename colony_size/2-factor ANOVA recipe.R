#script to run two-factor anova in survey package on univariate metrics from Swains


#load libraries
library(tidyverse)
library(survey)
library(car) #Anova w/ Type III SS function
library(emmeans) #post-hoc options for Anova 

rm(list=ls())

####Read in response variable dataset---------------
setwd("C:/github/swains/colony_size")
dat<-read.csv("MOSP_binned_tail_bins.csv", stringsAsFactors=FALSE) %>% mutate_if(is.character,as.factor) #modify to your response variable of interest
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

dat<-left_join(dat,NH, by = "DEPTH_BIN")

#Calculate survey weights (inverse proportion weighting)
w.df<-dat %>%
  group_by(OBS_YEAR,DEPTH_BIN,NH) %>%
  summarise(n = length(unique(SITE)))

w.df$sw<-w.df$NH/w.df$n #calculate survey weights for each site
dat.sw<-left_join(dat,w.df) #merge weights with site-level data
dat.sw$DEPTH_BIN <- as.factor(dat.sw$DEPTH_BIN)
head(dat.sw)
nrow(dat.sw) == n_distinct(dat.sw$SITE) #Check for duplicate sites

dat.sw<-filter(dat.sw, n!= "1") #remove strata that have less than 2 sites

dat.sw$Strat_conc<-as.factor(paste(dat.sw$OBS_YEAR, dat.sw$DEPTH_BIN,sep = "_"))
dat.des<-svydesign(id=~1, strata=~ Strat_conc, weights=~sw,data=dat.sw)




#Step 1: check whether we can use GLM or if non-parametric svyranktest is needed 

#test GLM assumptions of normality (shapiro.test) and homogeneity of variance (barlett's test) for each main effect (OBS_YEAR, DEPTH_BIN)
#replace "RV" with the name of your response variable column 
with(dat, tapply((RV), OBS_YEAR, shapiro.test)) # if p > 0.05 confirmed normal distribution among factor levels within OBS_YEAR
bartlett.test(RV ~ OBS_YEAR, dat) # p > 0.05 confirmed homogeneity of variance 

with(dat, tapply((RV), DEPTH_BIN, shapiro.test)) 
bartlett.test(RV ~ DEPTH_BIN, dat) 



#Step 2: if passed shapiro and bartlett test, good to proceed with GLM, with inverse-probability weighting and design-based standard errors (aka using survey package)
mod1<-svyglm(RV ~ DEPTH_BIN*OBS_YEAR, design=dat.des, family = gaussian()) #used gaussian for the size data (due to continuous, not count data)
Anova(mod1, type = 3, test.statistic = "F") 
emmeans(mod1, pairwise ~ OBS_YEAR) #post-hoc option for a significant main effect

plot(mod1) #check out residuals etc of the fitted model for any issues



#Step 2 alternate: if we didn't find parametric assumptions metin step 1, use this function instead
svyranktest(RV ~ DEPTH_BIN, design=dat.des, test=("KruskalWallis")) 
svyranktest(RV ~ OBS_YEAR, design=dat.des, test=("KruskalWallis"))
svyranktest(RV ~ DEPTH_BIN:OBS_YEAR, design=dat.des, test=("KruskalWallis")) 


#generate weighted means and design based standard errors for results resporting or making plots or your GLM model outputs
dat_weighted<-svyby(~RV,~OBS_YEAR,dat.des,svymean)