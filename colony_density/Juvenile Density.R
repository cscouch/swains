#Data do not include FUSP

# Using R version 4.1.0 (2021-05-18)

rm(list=ls())

library(dplyr)
library(tidyr)
library(survey)
library(multcomp)
library(emmeans)

#LOAD site-level data
site.data.gen2<-read.csv("T:/Benthic/Projects/Swains 2023 Benthic Analysis/Data/Swains_sitedata_GENUS.csv")
site.data.sp2<-read.csv("T:/Benthic/Projects/Swains 2023 Benthic Analysis/Data/Swains_sitedata_SPCODE.csv")
site.data.tax2<-write.csv("T:/Benthic/Projects/Swains 2023 Benthic Analysis/Data/Swains_sitedata_TAXONCODE.csv")

swa<-filter(site.data.gen2,new_MAX_DEPTH_M >=3) #subset sites less than 3m?

nrow(swa)

#read in sector-area file
sectors<-read.csv("T:/Benthic/Projects/Swains 2023 Benthic Analysis/Sectors-Strata-Areas.csv", stringsAsFactors=FALSE)
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
# site.sw$Strat_conc<-paste(site.sw$OBS_YEAR, site.sw$DB_RZ,sep = "_")
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


#Calculate regional mean and SE....of coral density?
# depth_mean<-svyby(~CORAL,~OBS_YEAR+DEPTH_BIN,des,svymean) # Dont have "CORAL" column
depth_mean<-svyby(~JuvColDen,~OBS_YEAR+DEPTH_BIN,des,svymean)
depth_mean_pocs<-svyby(~JuvColDen,~OBS_YEAR+DEPTH_BIN,des.pocs,svymean)
depth_mean_mosp<-svyby(~JuvColDen,~OBS_YEAR+DEPTH_BIN,des.mosp,svymean)
depth_mean_posp<-svyby(~JuvColDen,~OBS_YEAR+DEPTH_BIN,des.posp,svymean)
# depth_mean<-svyby(~JuvColDen,~OBS_YEAR+DEPTH_BIN,des,svymean)
depth_mean

#  Test fixed effects of region and year
modR<-svyglm(JuvColCount ~ OBS_YEAR*DEPTH_BIN, offset = TRANSECTAREA_ad, family = "poisson", design=des) # all significant
modR.pocs<-svyglm(JuvColCount ~ OBS_YEAR*DEPTH_BIN, offset = TRANSECTAREA_ad,family = "poisson", design=des.pocs) # Year and depth significant
modR.mosp<-svyglm(JuvColCount ~ OBS_YEAR*DEPTH_BIN, offset = TRANSECTAREA_ad,family = "poisson", design=des.mosp) # Depth and year:depth significant
modR.posp<-svyglm(JuvColCount ~ OBS_YEAR*DEPTH_BIN, offset = TRANSECTAREA_ad,family = "poisson", design=des.posp) # Depth significant

summary(modR); summary(modR.pocs); summary(modR.mosp); summary(modR.posp)
car::Anova(modR, type = 3, test.statistic = "F") #only year sig
car::Anova(modR.pocs, type = 3, test.statistic = "F") #year and the interaction sig
car::Anova(modR.mosp, type = 3, test.statistic = "F") #all factors sig
car::Anova(modR.posp, type = 3, test.statistic = "F") #

library(svydiags)
svyqqplot(JuvColDen~OBS_YEAR, design=des)
svystdres(modR,doplot=TRUE) 
bartlett.test(JuvColDen ~ OBS_YEAR, site.sw) #confirmed homogeneity of variance among factor levels within OBS_YEAR


# Post-hoc
# Year
y <- svyglm(JuvColCount ~ OBS_YEAR, offset = TRANSECTAREA_ad, family = "poisson", design = subset(des,DEPTH_BIN=="Shallow")) # All corals
emmeans(y, pairwise ~ OBS_YEAR) #post-hoc option for a significant main effect

y <- svyglm(JuvColCount ~ OBS_YEAR, offset = TRANSECTAREA_ad, family = "poisson", design = des.pocs) # POCS
summary(glht(y, mcp(OBS_YEAR="Tukey")))  # 2018 sig different than 2015 and 2023
emmeans(y, pairwise ~ OBS_YEAR) #post-hoc option for a significant main effect

y <- svyglm(JuvColCount ~ DEPTH_BIN, offset = TRANSECTAREA_ad, family = "poisson", design = des) # POCS
summary(glht(y, mcp(DEPTH_BIN="Tukey")))  # error
emmeans(y, pairwise ~ DEPTH_BIN) #post-hoc option for a significant main effect

summary(glht(modR.pocs, mcp(OBS_YEAR="Tukey"))) # 2023 sig different from 2015
summary(glht(modR.mosp, mcp(OBS_YEAR="Tukey"))) # 2023 and 2018 sig different than 2015
summary(glht(modR.posp, mcp(OBS_YEAR="Tukey"))) # 2023 and 2018 sig different than 2015


# Put all depth_mean dataframes together for plotting
depth_mean$Genus <- paste("Total")
depth_mean_pocs$Genus <- paste("Pocillopora")
depth_mean_mosp$Genus <- paste("Montipora")
depth_mean_posp$Genus <- paste("Porites")
depth_mean_tot <- rbind(depth_mean, depth_mean_pocs, depth_mean_mosp, depth_mean_posp)

# Plot Data
s <- ggplot(depth_mean_tot %>% filter(DEPTH_BIN == "Shallow"), 
            aes(x = OBS_YEAR, y = JuvColDen, group = Genus, fill = OBS_YEAR)) +
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = alpha(c("#440154FF","#22A884FF","#FDE725FF"))) +
  geom_errorbar(data = depth_mean_tot %>% filter(DEPTH_BIN == "Shallow"), 
                aes(ymin = JuvColDen-se, ymax = JuvColDen+se),
                width = .2) +
  facet_wrap(~Genus, nrow = 1) +
  guides(fill="none") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) +
  ggtitle("Shallow")

m <- ggplot(depth_mean_tot %>% filter(DEPTH_BIN == "Mid"), 
            aes(x = OBS_YEAR, y = JuvColDen, group = Genus, fill = OBS_YEAR)) +
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = alpha(c("#440154FF","#22A884FF","#FDE725FF"))) +
  geom_errorbar(data = depth_mean_tot %>% filter(DEPTH_BIN == "Mid"),
                aes(ymin = JuvColDen-se, ymax = JuvColDen+se),
                width = .1) +
  facet_wrap(~Genus, nrow = 1) +
  guides(fill="none") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())  +
  ggtitle("Mid")

d <- ggplot(depth_mean_tot %>% filter(DEPTH_BIN == "Deep"), 
            aes(x = OBS_YEAR, y = JuvColDen, group = Genus, fill = OBS_YEAR)) +
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = alpha(c("#440154FF","#22A884FF","#FDE725FF"))) +
  geom_errorbar(data = depth_mean_tot %>% filter(DEPTH_BIN == "Deep"),
                aes(ymin = JuvColDen-se, ymax = JuvColDen+se),
                width = .1) +
  facet_wrap(~Genus, nrow = 1) +
  guides(fill="none") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())  +
  ggtitle("Deep")

tot <- ggplot(depth_mean, 
              aes(x = OBS_YEAR, y = JuvColDen, fill = OBS_YEAR)) +
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = alpha(c("#440154FF","#22A884FF","#FDE725FF"))) +
  geom_errorbar(data = depth_mean,
                aes(ymin = JuvColDen-se, ymax = JuvColDen+se),
                width = .1) +
  facet_wrap(~DEPTH_BIN, nrow = 1) +
  guides(fill="none") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) 


install.packages("tidyselect")
ggarrange(s, m, d, 
          labels = c("Shallow", "Mid","Deep"),
          nrow = 3)

