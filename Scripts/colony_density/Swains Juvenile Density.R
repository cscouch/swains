# This script reads in summarized/prepped site-level coral juvenile data and runs svyglm models and plots for paper

# Using R version 4.1.0 (2021-05-18)

rm(list=ls())
dir = Sys.info()[7]
setwd(paste0("C:/Users/", dir, "/Documents/GitHub/swains/"))

library(tidyr)
library(dplyr)
library(survey)
library(multcomp)
library(emmeans)
library(ggplot2)
library(gridExtra)
library(sjstats)
library(pscl)
library(svydiags)

#Read in site-level data
site.data.gen2<-read.csv("Data/Swains_sitedata_GENUS.csv")


# FINAL DATA PREP ---------------------------------------------------------

#subset sites less than 3m because we did not survey 0-3m consistently across years and this habitat is very different than 3-6m
swa<-filter(site.data.gen2,new_MAX_DEPTH_M >=3) 
nrow(swa)

#read in sector-area file
sectors<-read.csv("Data/Sectors-Strata-Areas.csv", stringsAsFactors=FALSE)
swa_sa<-filter(sectors,ISLAND=="Swains")
swa_sa$DEPTH_BIN<-as.factor(swa_sa$DEPTH_BIN)

#Create a dataframe with the sum of NH (number of possible sites in a stratum) pooled across the sanctuary and open strata- we don't have enough sampling in open areas across time to separate into own strata
NH <- swa_sa %>%
  group_by(SEC_NAME, DEPTH_BIN)%>%
  #summarize(unique = NH)%>%
  group_by(DEPTH_BIN)%>%
  summarise(NH = sum(NH))

swa<-left_join(swa,NH) #merge juv data with new NH values pooled across the 2 swains sectors


#Calculate survey weights (inverse proportion weighting)
w.df<-swa %>%
  group_by(OBS_YEAR,DEPTH_BIN,NH) %>%
  summarise(n = length(unique(SITE)))

w.df$sw<-w.df$NH/w.df$n #calculate survey weights for each site

site.sw<-left_join(swa,w.df) #merge weights with site-level data
head(site.sw)


#Create contactenated Strata variable to specify as the design
site.sw$Strat_conc<-paste(site.sw$OBS_YEAR, site.sw$DEPTH_BIN,sep = "_")
site.sw$OBS_YEAR<-as.factor(site.sw$OBS_YEAR)

#Subset data for the 3 dominant taxa & total juveniles
site.sw.pocs <- site.sw %>% filter(GENUS_CODE == "POCS") # 
site.sw.mosp <- site.sw %>% filter(GENUS_CODE == "MOSP") # 
site.sw.posp <- site.sw %>% filter(GENUS_CODE == "POSP") # 
site.sw <- site.sw %>% filter(GENUS_CODE == "SSSS")

#remove strata that have less than 2 sites
site.sw.pocs <- site.sw.pocs %>% filter(n>1) # 
site.sw.mosp <- site.sw.mosp %>% filter(n>1) # 
site.sw.posp <- site.sw.posp %>% filter(n>1) # 
site.sw <- site.sw %>% filter(n>1)

summary(site.sw$n)
length(unique(site.sw$SITE)) #should be 50

#Check for duplicate sites- df should be empty
site.sw %>% 
  group_by(SITE) %>% 
  filter(n()>1)

head(site.sw)



# SURVEY MODELS -----------------------------------------------------------

#You should have a site-level dataframe with year, region, isalnd, sector, reef zone and depth bin, as well as a column(s) for your response variable and survey weights (sw)
#Now you are ready to use survey functions.

site.sw$OBS_YEAR<-as.factor(site.sw$OBS_YEAR)
site.sw$DEPTH_BIN<-as.factor(site.sw$DEPTH_BIN)

#Establish survey design with survey weights
des<-svydesign(id=~1, strata=~ Strat_conc, weights=~sw,data=site.sw)
des.pocs<-svydesign(id=~1, strata=~ Strat_conc, weights=~sw,data=site.sw.pocs)
des.mosp<-svydesign(id=~1, strata=~ Strat_conc, weights=~sw,data=site.sw.mosp)
des.posp<-svydesign(id=~1, strata=~ Strat_conc, weights=~sw,data=site.sw.posp)


#Use survey package to calculate mean SE and conduct statistical analyses
depth_mean<-svyby(~JuvColDen,~OBS_YEAR+DEPTH_BIN,des,svymean,na.rm.all=TRUE)
depth_mean_pocs<-svyby(~JuvColDen,~OBS_YEAR+DEPTH_BIN,des.pocs,svymean,na.rm.all=TRUE)
depth_mean_mosp<-svyby(~JuvColDen,~OBS_YEAR+DEPTH_BIN,des.mosp,svymean,na.rm.all=TRUE)
depth_mean_posp<-svyby(~JuvColDen,~OBS_YEAR+DEPTH_BIN,des.posp,svymean,na.rm.all=TRUE)
depth_mean


#Montipora 
#Shallow
site.sw.mosp.sh <- site.sw.mosp %>% filter(DEPTH_BIN == "Shallow") # 
des.mosp.sh<-svydesign(id=~1, strata=~ Strat_conc, weights=~sw,data=site.sw.mosp.sh)
modR.mosp.sh<-svyglm(JuvColCount ~ OBS_YEAR, offset = TRANSECTAREA_j, family = "poisson",design=des.mosp.sh) 
car::Anova(modR.mosp.sh, type = 3, test.statistic = "F") #nothing significant

#Mid
site.sw.mosp.mid <- site.sw.mosp %>% filter(DEPTH_BIN == "Mid") # 
des.mosp.mid<-svydesign(id=~1, strata=~ Strat_conc, weights=~sw,data=site.sw.mosp.mid)
modR.mosp.mid<-svyglm(JuvColCount ~ OBS_YEAR, offset = TRANSECTAREA_j, family = "poisson",design=des.mosp.mid) 
car::Anova(modR.mosp.mid, type = 3, test.statistic = "F") #nothing significant


#Deep
site.sw.mosp.deep <- site.sw.mosp %>% filter(DEPTH_BIN == "Deep") # 
des.mosp.deep<-svydesign(id=~1, strata=~ Strat_conc, weights=~sw,data=site.sw.mosp.deep)
modR.mosp.deep<-svyglm(JuvColCount ~ OBS_YEAR, offset = TRANSECTAREA_j, family = "poisson",design=des.mosp.deep) 
car::Anova(modR.mosp.deep, type = 3, test.statistic = "F") #nothing significant


#POCS
#Shallow
site.sw.pocs.sh <- site.sw.pocs %>% filter(DEPTH_BIN == "Shallow") # 
site.sw.pocs.sh$JuvColCount #lots of zeros

des.pocs.sh<-svydesign(id=~1, strata=~ Strat_conc, weights=~sw,data=site.sw.pocs.sh)
modR.pocs.sh<-svyglm(JuvColCount ~ OBS_YEAR, offset = TRANSECTAREA_j, family = "poisson",design=des.pocs.sh) 
car::Anova(modR.pocs.sh, type = 3, test.statistic = "F") #nothing significant

#Mid
site.sw.pocs.mid <- site.sw.pocs %>% filter(DEPTH_BIN == "Mid") # 
site.sw.pocs.mid$JuvColCount #lots of zeros
des.pocs.mid<-svydesign(id=~1, strata=~ Strat_conc, weights=~sw,data=site.sw.pocs.mid)
modR.pocs.mid<-svyglm(JuvColCount ~ OBS_YEAR, offset = TRANSECTAREA_j, family = "poisson",design=des.pocs.mid) 
car::Anova(modR.pocs.mid, type = 3, test.statistic = "F") #nothing significant
AIC(modR.mosp.mid)

#Try out tweedie model given strong zero inflation
mod.tw = svyglm(JuvColCount ~ OBS_YEAR, offset = TRANSECTAREA_j, design = des.pocs.mid,
                family = statmod::tweedie(var.power = 1.2, link.power = 0), maxit = 500)
AIC(mod.tw) #tweedie model has lower AIC than quasipoisson- moving forward with tweedie

car::Anova(mod.tw, type = 3, test.statistic = "F") #year is significant
emmeans(mod.tw, specs = pairwise~OBS_YEAR, adjust = "none") #posthoc tests - add letters to figure


#Deep-only 1 site had POCS juveniles- do not proceed with models
site.sw.pocs.deep <- site.sw.pocs %>% filter(DEPTH_BIN == "Deep") # 
site.sw.pocs.deep$JuvColCount #lots of zeros



#POSP
#Shallow
bartlett.test(JuvColDen ~ Strat_conc, subset(site.sw.posp, DEPTH_BIN == "Shallow")) 
site.sw.posp.sh <- site.sw.posp %>% filter(DEPTH_BIN == "Shallow") # 
site.sw.posp.sh$JuvColCount #lots of zeros

des.posp.sh<-svydesign(id=~1, strata=~ Strat_conc, weights=~sw,data=site.sw.posp.sh)
modR.posp.sh<-svyglm(JuvColCount ~ OBS_YEAR, offset = TRANSECTAREA_j, family = "quasipoisson",design=des.posp.sh) 
car::Anova(modR.posp.sh, type = 3, test.statistic = "F") #nothing significant
AIC(modR.posp.sh)

#Try out tweedie model given strong zero inflation
mod.tw = svyglm(JuvColCount ~ OBS_YEAR, offset = TRANSECTAREA_j, design = des.posp.sh,
                family = statmod::tweedie(var.power = 1.2, link.power = 0), maxit = 500)
AIC(mod.tw) #tweedie model has lower AIC than quasipoisson- moving forward with tweedie

car::Anova(mod.tw, type = 3, test.statistic = "F") #year is significant
emmeans(mod.tw, specs = pairwise~OBS_YEAR, adjust = "none") #posthoc tests - add letters to figure


#Mid - can't run test (only 1 site has value > 0)
bartlett.test(sqrt(JuvColDen) ~ Strat_conc, subset(site.sw.posp, DEPTH_BIN == "Mid")) 
site.sw.posp.mid <- site.sw.posp %>% filter(DEPTH_BIN == "Mid") # 
site.sw.posp.mid$JuvColCount #only one site has value >0


#Deep -too few sites with too many zeros- do not run models
bartlett.test(JuvColDen ~ Strat_conc, subset(site.sw.posp, DEPTH_BIN == "Deep")) 
site.sw.posp.deep <- site.sw.posp %>% filter(DEPTH_BIN == "Deep") # 
site.sw.posp.deep$JuvColCount #only one site has value >0



# PLOTTING ----------------------------------------------------------------
#Plot 3 dominant taxa for all depth bins 

# Put all depth_mean dataframes together for plotting
depth_mean_pocs$Genus <- paste("Pocillopora")
depth_mean_mosp$Genus <- paste("Montipora")
depth_mean_posp$Genus <- paste("Porites")
depth_mean_tot <- rbind(depth_mean_pocs, depth_mean_mosp, depth_mean_posp)

#Shallow plot
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
  scale_y_continuous(expand = c(0,0), limits = c(-0.1,6.5))+
  geom_text(aes(x=OBS_YEAR,y=JuvColDen+se,label=c("","","","","","","a","b","b")),
            position = position_dodge(),
            vjust = -0.5)

#Mid plot
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


#Deep plot
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


ytitle= "Mean Juvenile Colonies m-2"
xtitle= "Year"

juv.plot<-grid.arrange(s+ggtitle("Shallow"), m + ggtitle("Mid"), d+ ggtitle("Deep"), 
          left = ytitle,
          bottom = xtitle,
          nrow = 3)

#ggsave(plot=juv.plot,file="T:/Benthic/Projects/Swains 2023 Benthic Analysis/Plots/FigureS3.jpg",width=10,height=8)

ggsave(plot=juv.plot,file="FigureS3.jpg",width=10,height=8)
