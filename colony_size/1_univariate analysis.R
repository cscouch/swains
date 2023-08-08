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
####Read in binned data---------------
setwd("C:/github/swains/colony_size")
MOSP<-read.csv("MOSP_binned_tail_bins.csv", stringsAsFactors=FALSE) %>% mutate_if(is.character,as.factor)
POCS<-read.csv("POCS_binned_tail_bins.csv", stringsAsFactors=FALSE) %>% mutate_if(is.character,as.factor)
POSP<-read.csv("POSP_binned_tail_bins.csv", stringsAsFactors=FALSE) %>% mutate_if(is.character,as.factor)

MOSP$OBS_YEAR <- as.factor(MOSP$OBS_YEAR)
POCS$OBS_YEAR <- as.factor(POCS$OBS_YEAR)
POSP$OBS_YEAR <- as.factor(POSP$OBS_YEAR)



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
MOSP.des<-svydesign(id=~1, strata=~ Strat_conc, weights=~sw,data=MOSP.sw)


#POSP
w.df<-POCS %>%
  group_by(OBS_YEAR,REEF_ZONE,DEPTH_BIN,NH) %>%
  summarise(n = length(unique(SITE)))

w.df$sw<-w.df$NH/w.df$n #calculate survey weights for each site
POCS.sw<-left_join(POCS,w.df) #merge weights with site-level data
POCS.sw$DEPTH_BIN <- as.factor(POCS.sw$DEPTH_BIN)
head(POCS.sw)
nrow(POCS.sw) == n_distinct(POCS.sw$SITE) #Check for duplicate sites

POCS.sw$Strat_conc<-as.factor(paste(POCS.sw$OBS_YEAR, POCS.sw$DEPTH_BIN,sep = "_"))
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

POSP.sw$Strat_conc<-as.factor(paste(POSP.sw$OBS_YEAR, POSP.sw$DEPTH_BIN,sep = "_"))
POSP.des<-svydesign(id=~1, strata=~ Strat_conc, weights=~sw,data=POSP.sw)




#### MEAN SIZE ANALYSIS: two-factor (YEAR *DEPTH) using site-level data-------------------------------

#Test fixed effects of region and year
#MOSP
MOSP1<-svyglm(QJuv.R ~ DEPTH_BIN*OBS_YEAR, design=MOSP.des)
anova(MOSP1) 
#NS

MOSP2<-svyglm(Q10.R ~ DEPTH_BIN*OBS_YEAR, design=MOSP.des)
anova(MOSP2)
#NS

MOSP3<-svyglm(QMed.R ~ DEPTH_BIN*OBS_YEAR, design=MOSP.des)
anova(MOSP3)
#NS

MOSP4<-svyglm(Q90.R ~ DEPTH_BIN*OBS_YEAR, design=MOSP.des)
anova(MOSP4) #depth bin is significant
emmeans(MOSP4, pairwise ~ DEPTH_BIN) #shallow differs from mod for the large size class only



#POCS
POCS1<-svyglm(QJuv.R ~ DEPTH_BIN*OBS_YEAR, design=POCS.des)
anova(POCS1) 
#NS

POCS2<-svyglm(Q10.R ~ DEPTH_BIN*OBS_YEAR, design=POCS.des)
anova(POCS2) #OBS year is barely significant
emmeans(POCS2, pairwise ~ OBS_YEAR) #but no differences in post-hocs with adjusted p-values

POCS3<-svyglm(QMed.R ~ DEPTH_BIN*OBS_YEAR, design=POCS.des)
anova(POCS3)
#NS

POCS4<-svyglm(Q90.R ~ DEPTH_BIN*OBS_YEAR, design=POCS.des)
anova(POCS4) #depth bin is significant
emmeans(POCS4, pairwise ~ DEPTH_BIN) #shallow differs from mod for the large size class only




#POSP
POSP1<-svyglm(QJuv.R ~ DEPTH_BIN*OBS_YEAR, design=POSP.des)
anova(POSP1) 
#NS

POSP2<-svyglm(Q10.R ~ DEPTH_BIN*OBS_YEAR, design=POSP.des)
anova(POSP2) 
#NS

POSP3<-svyglm(QMed.R ~ DEPTH_BIN*OBS_YEAR, design=POSP.des)
anova(POSP3)
#NS

POSP4<-svyglm(Q90.R ~ DEPTH_BIN*OBS_YEAR, design=POSP.des)
anova(POSP4) #OBS_YEAR is significant
emmeans(POSP4, pairwise ~ OBS_YEAR) #hmmm.....seems misleading...no data?




####GRAPHIING -------------------------------
setwd("C:/github/swains/colony_size")

#Calculate means by depth as this was the only significant factor
#MOSP
M1<-svyby(~QJuv.R,~DEPTH_BIN,MOSP.des,svymean)
M2<-svyby(~Q10.R,~DEPTH_BIN,MOSP.des,svymean)
M3<-svyby(~QMed.R,~DEPTH_BIN,MOSP.des,svymean)
M4<-svyby(~Q90.R,~DEPTH_BIN,MOSP.des,svymean)
names(M1)<-c("Depth_Bin", "mean", "se")
names(M2)<-c("Depth_Bin", "mean", "se")
names(M3)<-c("Depth_Bin", "mean", "se")
names(M4)<-c("Depth_Bin", "mean", "se")
M1$BIN <- as.factor("Juv")
M2$BIN <- as.factor("Sm")
M3$BIN <- as.factor("Med")
M4$BIN <- as.factor("Lg")
MOSP_summary <- rbind(M1,M2,M3,M4)

MOSP_summary$group <- c( "", "", "",  
                     "", "", "",
                     "", "", "",
                     "AB", "B", "A")

#bar plot of mean size
g1 <- ggplot(MOSP_summary, aes(x = BIN, y = mean, fill = Depth_Bin)) + 
  geom_bar(position="dodge", stat="identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2,position=position_dodge(.9)) +
  geom_text(data=MOSP_summary,aes(x = BIN, y = mean+se, label=group), size = 4, position=position_dodge(0.9),  vjust = -0.65) +
  scale_fill_viridis_d() +
  theme_ipsum() +
  ylab("Proportion")+
  xlab("Size Bin")+
  ggtitle("MOSP") #name to appropriate taxa

g1
ggsave ("MOSP_final.jpeg", width =5, height =5, units = 'in')


#POCS ---summary for two factor crossed
P1<-svyby(~QJuv.R,~OBS_YEAR+DEPTH_BIN,POCS.des,svymean)
P2<-svyby(~Q10.R,~OBS_YEAR+DEPTH_BIN,POCS.des,svymean)
P3<-svyby(~QMed.R,~OBS_YEAR+DEPTH_BIN,POCS.des,svymean)
P4<-svyby(~Q90.R,~OBS_YEAR+DEPTH_BIN,POCS.des,svymean)
names(P1)<-c("Year", "Depth_Bin", "mean", "se")
names(P2)<-c("Year", "Depth_Bin", "mean", "se")
names(P3)<-c("Year", "Depth_Bin", "mean", "se")
names(P4)<-c("Year", "Depth_Bin", "mean", "se")
P1$BIN <- as.factor("Juv")
P2$BIN <- as.factor("Sm")
P3$BIN <- as.factor("Med")
P4$BIN <- as.factor("Lg")
POCS_summary <- rbind(P1,P2,P3,P4)


g2 <- ggplot(POCS_summary, aes(x = BIN, y = mean, fill = Year)) + 
  geom_bar(position="dodge", stat="identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2,position=position_dodge(.9)) +
  facet_wrap("Depth_Bin",  ncol =1, scales = "fixed")+
  scale_fill_viridis_d() +
  theme_ipsum() +
  ylab("Proportion")+
  xlab("Size Bin")+
  ggtitle("POCS") #name to appropriate taxa
g2


#POCS ---summary for depth only
M1<-svyby(~QJuv.R,~DEPTH_BIN,POCS.des,svymean)
M2<-svyby(~Q10.R,~DEPTH_BIN,POCS.des,svymean)
M3<-svyby(~QMed.R,~DEPTH_BIN,POCS.des,svymean)
M4<-svyby(~Q90.R,~DEPTH_BIN,POCS.des,svymean)
names(M1)<-c("Depth_Bin", "mean", "se")
names(M2)<-c("Depth_Bin", "mean", "se")
names(M3)<-c("Depth_Bin", "mean", "se")
names(M4)<-c("Depth_Bin", "mean", "se")
M1$BIN <- as.factor("Juv")
M2$BIN <- as.factor("Sm")
M3$BIN <- as.factor("Med")
M4$BIN <- as.factor("Lg")
POCS_summary <- rbind(M1,M2,M3,M4)

POCS_summary$group <- c( "", "", "",  
                         "", "", "",
                         "", "", "",
                         "AB", "B", "A")

#bar plot of mean size
g3 <- ggplot(POCS_summary, aes(x = BIN, y = mean, fill = Depth_Bin)) + 
  geom_bar(position="dodge", stat="identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2,position=position_dodge(.9)) +
  geom_text(data=POCS_summary,aes(x = BIN, y = mean+se, label=group), size = 4, position=position_dodge(0.9),  vjust = -0.65) +
  scale_fill_viridis_d() +
  theme_ipsum() +
  ylab("Proportion")+
  xlab("Size Bin")+
  ggtitle("POCS") #name to appropriate taxa

g3
ggsave ("POCS_final.jpeg", width =5, height =5, units = 'in')


#POSP--summary for obs year only
P1<-svyby(~QJuv.R,~OBS_YEAR,POSP.des,svymean)
P2<-svyby(~Q10.R,~OBS_YEAR,POSP.des,svymean)
P3<-svyby(~QMed.R,~OBS_YEAR,POSP.des,svymean)
P4<-svyby(~Q90.R,~OBS_YEAR,POSP.des,svymean)
names(P1)<-c("Year", "mean", "se")
names(P2)<-c("Year", "mean", "se")
names(P3)<-c("Year", "mean", "se")
names(P4)<-c("Year", "mean", "se")
P1$BIN <- as.factor("Juv")
P2$BIN <- as.factor("Sm")
P3$BIN <- as.factor("Med")
P4$BIN <- as.factor("Lg")
POSP_summary <- rbind(P1,P2,P3,P4)

#bar plot
g4 <- ggplot(POSP_summary, aes(x = BIN, y = mean, fill = Year)) + 
  geom_bar(position="dodge", stat="identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2,position=position_dodge(.9)) +
  #geom_text(data=POCS_summary,aes(x = BIN, y = mean+se, label=group), size = 4, position=position_dodge(0.9),  vjust = -0.65) +
  scale_fill_viridis_d() +
  theme_ipsum() +
  ylab("Proportion")+
  xlab("Size Bin")+
  ggtitle("POSP") #name to appropriate taxa

g4
ggsave ("POSP_final.jpeg", width =5, height =5, units = 'in')





