#Data do not include FUSP

# Using R version 4.1.0 (2021-05-18)

# This script will clean the raw colony-level NOAA NCRMP data
rm(list=ls())
dir = Sys.info()[7]
setwd(paste0("C:/Users/", dir, "/Documents/GitHub/swains/"))

library(dplyr)
library(tidyr)
library(survey)
library(multcomp)
library(emmeans)
library(ggplot2)
library(gridExtra)
library(sjstats)
library(pscl)
library(svydiags)

#LOAD site-level data
site.data.gen2<-read.csv("Data/Swains_sitedata_GENUS.csv")

swa<-filter(site.data.gen2,new_MAX_DEPTH_M >=3) #subset sites less than 3m?

nrow(swa)

#read in sector-area file
sectors<-read.csv("Data/Sectors-Strata-Areas.csv", stringsAsFactors=FALSE)
swa_sa<-filter(sectors,ISLAND=="Swains")
swa_sa$DEPTH_BIN<-as.factor(swa_sa$DEPTH_BIN)

NH <- swa_sa %>%
  group_by(SEC_NAME, DEPTH_BIN)%>%
  summarize(unique = NH)%>%
  group_by(DEPTH_BIN)%>%
  summarise(NH = sum(unique))

swa<-left_join(swa,NH) #merge demo data with new NH values pooled across the 2 swains sectors



#Calculate survey weights (inverse proportion weighting)
w.df<-swa %>%
  group_by(OBS_YEAR,DEPTH_BIN,NH) %>%
  summarise(n = length(unique(SITE)))

w.df$sw<-w.df$NH/w.df$n #calculate survey weights for each site

site.sw<-left_join(swa,w.df) #merge weights with site-level data
head(site.sw)

nrow(site.sw) #should be 50


#Create contactenated Strata variable 
site.sw$Strat_conc<-paste(site.sw$OBS_YEAR, site.sw$DEPTH_BIN,sep = "_")
site.sw$OBS_YEAR<-as.factor(site.sw$OBS_YEAR)

#Subset data
site.sw.pocs <- site.sw %>% filter(GENUS_CODE == "POCS") # above code isn't 50 so here I'm looking at just total coral
site.sw.mosp <- site.sw %>% filter(GENUS_CODE == "MOSP") # 
site.sw.posp <- site.sw %>% filter(GENUS_CODE == "POSP") # 
site.sw <- site.sw %>% filter(GENUS_CODE == "SSSS")

#remove strata that have less than 2 sites
site.sw<-subset(site.sw,n>1)
summary(site.sw$n)
length(unique(site.sw$SITE))

site.sw.pocs<-subset(site.sw.pocs,n>1)
site.sw.mosp<-subset(site.sw.mosp,n>1)
site.sw.posp<-subset(site.sw.posp,n>1)

#Check for duplicate sites- df should be empty
site.sw %>% 
  group_by(SITE) %>% 
  filter(n()>1)

head(site.sw)

#You should have a site-level dataframe with year, region, isalnd, sector, reef zone and depth bin, as well as a column(s) for your response variable and survey weights (sw)
#Now you are ready to use survey functions.

#Use survey package to calculate mean SE and conduct statistical analyses
site.sw$OBS_YEAR<-as.factor(site.sw$OBS_YEAR)
site.sw$DEPTH_BIN<-as.factor(site.sw$DEPTH_BIN)

#Establish survey design with survey weights
des<-svydesign(id=~1, strata=~ Strat_conc, weights=~sw,data=site.sw)
des.pocs<-svydesign(id=~1, strata=~ Strat_conc, weights=~sw,data=site.sw.pocs)
des.mosp<-svydesign(id=~1, strata=~ Strat_conc, weights=~sw,data=site.sw.mosp)
des.posp<-svydesign(id=~1, strata=~ Strat_conc, weights=~sw,data=site.sw.posp)


# Analyze juvenile density 1-5cm ------------------------------------------
depth_mean<-svyby(~JuvColDen,~OBS_YEAR+DEPTH_BIN,des,svymean,na.rm.all=TRUE)
depth_mean_pocs<-svyby(~JuvColDen,~OBS_YEAR+DEPTH_BIN,des.pocs,svymean,na.rm.all=TRUE)
depth_mean_mosp<-svyby(~JuvColDen,~OBS_YEAR+DEPTH_BIN,des.mosp,svymean,na.rm.all=TRUE)
depth_mean_posp<-svyby(~JuvColDen,~OBS_YEAR+DEPTH_BIN,des.posp,svymean,na.rm.all=TRUE)
# depth_mean<-svyby(~JuvColDen,~OBS_YEAR+DEPTH_BIN,des,svymean)
depth_mean

#Total 
#Shallow
with(subset(site.sw, DEPTH_BIN == "Shallow"), tapply((JuvColDen), OBS_YEAR, shapiro.test))
bartlett.test(JuvColDen ~ Strat_conc, subset(site.sw, DEPTH_BIN == "Shallow")) 
site.sw.sh <- site.sw %>% filter(DEPTH_BIN == "Shallow") # 
des<-svydesign(id=~1, strata=~ Strat_conc, weights=~sw,data=site.sw.sh)

modR<-svyglm(JuvColDen ~ OBS_YEAR,  design=des) 
car::Anova(modR, type = 3, test.statistic = "F") #nothing significant



#MOSP
#Shallow
with(subset(site.sw.mosp, DEPTH_BIN == "Shallow"), tapply((JuvColDen), OBS_YEAR, shapiro.test))
bartlett.test(JuvColDen ~ Strat_conc, subset(site.sw.mosp, DEPTH_BIN == "Shallow")) 
modR.mosp<-svyglm(JuvColDen ~ OBS_YEAR,  design=subset(des.mosp,DEPTH_BIN=="Shallow")) 
car::Anova(modR.mosp, type = 3, test.statistic = "F") #nothing significant

#Mid
with(subset(site.sw.mosp, DEPTH_BIN == "Mid"), tapply((JuvColDen), OBS_YEAR, shapiro.test))
bartlett.test(JuvColDen ~ Strat_conc, subset(site.sw.mosp, DEPTH_BIN == "Mid")) 
modR.mosp<-svyglm(JuvColDen ~ OBS_YEAR, design=subset(des.mosp,DEPTH_BIN=="Mid")) # Depth and year:depth significant
car::Anova(modR.mosp, type = 3, test.statistic = "F") #nothing significant


#Deep
with(subset(site.sw.mosp, DEPTH_BIN == "Deep"), tapply((JuvColDen), OBS_YEAR, shapiro.test))
bartlett.test(JuvColDen ~ Strat_conc, subset(site.sw.mosp, DEPTH_BIN == "Deep")) 
modR.mosp<-svyglm(JuvColDen ~ OBS_YEAR, design=subset(des.mosp,DEPTH_BIN=="Deep")) 
car::Anova(modR.mosp, type = 3, test.statistic = "F") #nothing significant


#POCS
#Shallow
with(subset(site.sw.pocs, DEPTH_BIN == "Shallow"), tapply((sqrt(JuvColDen)), OBS_YEAR, shapiro.test))
bartlett.test(sqrt(JuvColDen) ~ Strat_conc, subset(site.sw.pocs, DEPTH_BIN == "Shallow")) 
site.sw.pocs.sh <- site.sw.pocs %>% filter(DEPTH_BIN == "Shallow") # 
des.pocs.sh<-svydesign(id=~1, strata=~ Strat_conc, weights=~sw,data=site.sw.pocs.sh)
modR.pocs<-svyglm(sqrt(JuvColDen) ~ OBS_YEAR,  design=des.pocs.sh) 
car::Anova(modR.pocs, type = 3, test.statistic = "F") #nothing significant

#Mid
with(subset(site.sw.mosp, DEPTH_BIN == "Mid"), tapply((JuvColDen), OBS_YEAR, shapiro.test))
bartlett.test(sqrt(JuvColDen) ~ Strat_conc, subset(site.sw.pocs, DEPTH_BIN == "Mid")) 
site.sw.pocs.mid <- site.sw.pocs %>% filter(DEPTH_BIN == "Mid") # 
des.pocs.mid<-svydesign(id=~1, strata=~ Strat_conc, weights=~sw,data=site.sw.pocs.mid)
modR.pocs<-svyglm(JuvColDen ~ OBS_YEAR, design=subset(des.mosp,DEPTH_BIN=="Mid")) 
car::Anova(modR.pocs, type = 3, test.statistic = "F") #nothing significant


#Deep- use non parametric stats
with(subset(site.sw.pocs, DEPTH_BIN == "Deep"), tapply((JuvColDen), OBS_YEAR, shapiro.test))
bartlett.test(JuvColDen ~ Strat_conc, subset(site.sw.pocs, DEPTH_BIN == "Deep")) 
#svyranktest(JuvColDen ~ OBS_YEAR, design=subset(des.pocs,DEPTH_BIN=="Deep"), test=("KruskalWallis")) #not enough colonies to run model (only 1 site across all years had POCS juves)

#Take home: sqrt transform SSSS, MOSP, POCS, but POSP use non-parametric stats


#POSP
#Shallow
with(subset(site.sw.posp, DEPTH_BIN == "Shallow"), tapply((sqrt(JuvColDen)), OBS_YEAR, shapiro.test))
bartlett.test(sqrt(JuvColDen) ~ Strat_conc, subset(site.sw.posp, DEPTH_BIN == "Shallow")) 

#
site.sw.posp.sh <- site.sw.pocs %>% filter(DEPTH_BIN == "Shallow") # 
des.posp.sh<-svydesign(id=~1, strata=~ Strat_conc, weights=~sw,data=site.sw.posp.sh)

svyranktest(JuvColDen ~ OBS_YEAR, design=des.posp.sh, test=("KruskalWallis"))
# svyranktest(JuvColDen ~ OBS_YEAR, design=subset(des.posp,DEPTH_BIN=="Mid"), test=("KruskalWallis"))#not enough colonies to run model (only 1 site across all years had POSP juves)
svyranktest(JuvColDen ~ OBS_YEAR, design=subset(des.posp,DEPTH_BIN=="Deep"), test=("KruskalWallis"))

dunnTest(JuvColDen~OBS_YEAR,data=subset(site.sw.posp,DEPTH_BIN=="Shallow"))


l<-c("2015","2018","2023")
ps<-matrix(NA,3,3)
dimnames(ps)<-list(l,l)
for(i in 1:2){
  for(j in (i+1):3){
    ps[i,j]<-svyranktest(JuvColDen~OBS_YEAR, subset(des.posp, OBS_YEAR %in% l[c(i,j)]))$p.value
  }
}
ps

p.adjust(ps[!is.na(ps)],method="hochberg")


# Put all depth_mean dataframes together for plotting
#depth_mean$Genus <- paste("Total")
depth_mean_pocs$Genus <- paste("Pocillopora")
depth_mean_mosp$Genus <- paste("Montipora")
depth_mean_posp$Genus <- paste("Porites")
depth_mean_tot <- rbind(depth_mean_pocs, depth_mean_mosp, depth_mean_posp)

# Plot Data
# depth_mean_tot$sig<-c("","ns","","","ns","","a","b","b",
#                       "","ns","","","ns","","a","a","b",
#                       "","ns","","a","b","b","a","b","a")

s <- ggplot(depth_mean_tot %>% filter(DEPTH_BIN == "Shallow"),
            aes(x = OBS_YEAR, y = JuvColDen, group = Genus, fill = OBS_YEAR)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = alpha(c("#FF8C00","#44BB99","#AA4499"))) +
  geom_errorbar(data = depth_mean_tot %>% filter(DEPTH_BIN == "Shallow"),
                aes(ymin = JuvColDen-se, ymax = JuvColDen+se),
                width = .2) +
  facet_wrap(~Genus, nrow = 1) +
  guides(fill="none") +
  
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())+
  scale_y_continuous(expand = c(0,0), limits = c(-0.1,6.5)) 
  # geom_text(aes(x=OBS_YEAR,y=JuvColDen+se,label=c("","ns","","","ns","","a","b","b")),
  #           position = position_dodge(),
  #           vjust = -0.5) 

m <- ggplot(depth_mean_tot %>% filter(DEPTH_BIN == "Mid"), 
            aes(x = OBS_YEAR, y = JuvColDen, group = Genus, fill = OBS_YEAR)) +
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = alpha(c("#FF8C00","#44BB99","#AA4499"))) +
  geom_errorbar(data = depth_mean_tot %>% filter(DEPTH_BIN == "Mid"),
                aes(ymin = JuvColDen-se, ymax = JuvColDen+se),
                width = .1) +
  facet_wrap(~Genus, nrow = 1) +
  guides(fill="none") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())+
  scale_y_continuous(expand = c(0,0), limits = c(-0.1,6.5))
  # geom_text(aes(x=OBS_YEAR,y=JuvColDen+se,label=c("","ns","","","ns","","a","a","b")),
  #           position = position_dodge(),
  #           vjust = -0.5) 

d <- ggplot(depth_mean_tot %>% filter(DEPTH_BIN == "Deep"), 
            aes(x = OBS_YEAR, y = JuvColDen, group = Genus, fill = OBS_YEAR)) +
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = alpha(c("#FF8C00","#44BB99","#AA4499"))) +
  geom_errorbar(data = depth_mean_tot %>% filter(DEPTH_BIN == "Deep"),
                aes(ymin = JuvColDen-se, ymax = JuvColDen+se),
                width = .1) +
  facet_wrap(~Genus, nrow = 1) +
  guides(fill="none") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())+
  scale_y_continuous(expand = c(0,0), limits = c(-0.1,6.5))
  # geom_text(aes(x=OBS_YEAR,y=JuvColDen+se,label=c("a","b","b","","ns","","a","b","a")),
  #           position = position_dodge(),
  #           vjust = -0.5)   

tot <- ggplot(depth_mean, 
              aes(x = OBS_YEAR, y = JuvColDen, fill = OBS_YEAR)) +
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = alpha(c("#FF8C00","#44BB99","#AA4499"))) +
  geom_errorbar(data = depth_mean,
                aes(ymin = JuvColDen-se, ymax = JuvColDen+se),
                width = .1) +
  facet_wrap(~DEPTH_BIN, nrow = 1) +
  guides(fill="none") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) 

ytitle= "Mean Juvenile Colonies m-2"
xtitle= "Year"

juv.plot<-grid.arrange(s+ggtitle("Shallow"), m + ggtitle("Mid"), d+ ggtitle("Deep"), 
          left = ytitle,
          bottom = xtitle,
          nrow = 3)

ggsave(plot=juv.plot,file="T:/Benthic/Projects/Swains 2023 Benthic Analysis/Plots/FigureS3.jpg",width=10,height=8)

