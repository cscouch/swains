
# Using R version 4.1.0 (2021-05-18)

rm(list=ls())

library(dplyr)
library(tidyselect)
library(ggpubr)
library(survey)
library(emmeans)

#LOAD site-level data
site.data.gen2<-read.csv("T:/Benthic/Projects/Swains 2023 Benthic Analysis/Data/Swains_sitedata_GENUS.csv")
site.data.sp2<-read.csv("T:/Benthic/Projects/Swains 2023 Benthic Analysis/Data/Swains_sitedata_SPCODE.csv")
site.data.tax2<-read.csv("T:/Benthic/Projects/Swains 2023 Benthic Analysis/Data/Swains_sitedata_TAXONCODE.csv")
swan <- read.csv("T:/Benthic/Projects/Swains 2023 Benthic Analysis/Data/Swains_sitedata_TAXONCODE_MORPH.csv")
swam <- read.csv("T:/Benthic/Projects/Swains 2023 Benthic Analysis/Data/Swains_sitedata_TAXONCODE_MORPH.csv")
swa<-filter(site.data.gen2,new_MAX_DEPTH_M >=3) #subset sites less than 3m?

nrow(swa)

#read in sector-area file
sectors<-read.csv("T:/Benthic/Projects/Swains 2023 Benthic Analysis/Sectors-Strata-Areas.csv", stringsAsFactors=FALSE)
swa_sa<-filter(sectors,ISLAND=="Swains")

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


#### Genus level adult density ####
#Subset data
site.sw.pocs <- site.sw %>% filter(GENUS_CODE == "POCS") # above code isn't 50 so here I'm looking at just total coral
site.sw.mosp <- site.sw %>% filter(GENUS_CODE == "MOSP") # 
site.sw.posp <- site.sw %>% filter(GENUS_CODE == "POSP") # 
site.sw <- site.sw %>% filter(GENUS_CODE == "SSSS")

# #Subset data by species and depth 
# site.sw.pocs.s <- site.sw %>% filter(GENUS_CODE == "POCS" & DEPTH_BIN == "Shallow") # above code isn't 50 so here I'm looking at just total coral
# site.sw.mosp.s <- site.sw %>% filter(GENUS_CODE == "MOSP" & DEPTH_BIN == "Shallow") # 
# site.sw.posp.s <- site.sw %>% filter(GENUS_CODE == "POSP" & DEPTH_BIN == "Shallow") # 
# site.sw.s <- site.sw %>% filter(GENUS_CODE == "SSSS" & DEPTH_BIN == "Shallow")
# 
# site.sw.pocs.m <- site.sw %>% filter(GENUS_CODE == "POCS" & DEPTH_BIN == "Mid") # above code isn't 50 so here I'm looking at just total coral
# site.sw.mosp.m <- site.sw %>% filter(GENUS_CODE == "MOSP" & DEPTH_BIN == "Mid") # 
# site.sw.posp.m <- site.sw %>% filter(GENUS_CODE == "POSP" & DEPTH_BIN == "Mid") # 
# site.sw.m <- site.sw %>% filter(GENUS_CODE == "SSSS" & DEPTH_BIN == "Mid")
# 
# site.sw.pocs.d <- site.sw %>% filter(GENUS_CODE == "POCS" & DEPTH_BIN == "Deep") # above code isn't 50 so here I'm looking at just total coral
# site.sw.mosp.d <- site.sw %>% filter(GENUS_CODE == "MOSP" & DEPTH_BIN == "Deep") # 
# site.sw.posp.d <- site.sw %>% filter(GENUS_CODE == "POSP" & DEPTH_BIN == "Deep") # 
# site.sw.d <- site.sw %>% filter(GENUS_CODE == "SSSS" & DEPTH_BIN == "Deep")

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
site.sw.pocs$OBS_YEAR<-as.factor(site.sw.pocs$OBS_YEAR)
site.sw.mosp$OBS_YEAR<-as.factor(site.sw.mosp$OBS_YEAR)
site.sw.posp$OBS_YEAR<-as.factor(site.sw.posp$OBS_YEAR)

#Establish survey design with survey weights
des<-svydesign(id=~1, strata=~ Strat_conc, weights=~sw,data=site.sw)
des.pocs<-svydesign(id=~1, strata=~ Strat_conc, weights=~sw,data=site.sw.pocs)
des.mosp<-svydesign(id=~1, strata=~ Strat_conc, weights=~sw,data=site.sw.mosp)
des.posp<-svydesign(id=~1, strata=~ Strat_conc, weights=~sw,data=site.sw.posp)


# #Establish survey design with survey weights including depth
# des.s<-svydesign(id=~1, strata=~ Strat_conc, weights=~sw,data=site.sw.s)
# des.pocs.s<-svydesign(id=~1, strata=~ Strat_conc, weights=~sw,data=site.sw.pocs.s)
# des.mosp.s<-svydesign(id=~1, strata=~ Strat_conc, weights=~sw,data=site.sw.mosp.s)
# des.posp.s<-svydesign(id=~1, strata=~ Strat_conc, weights=~sw,data=site.sw.posp.s)
# 
# des.m<-svydesign(id=~1, strata=~ Strat_conc, weights=~sw,data=site.sw.m)
# des.pocs.m<-svydesign(id=~1, strata=~ Strat_conc, weights=~sw,data=site.sw.pocs.m)
# des.mosp.m<-svydesign(id=~1, strata=~ Strat_conc, weights=~sw,data=site.sw.mosp.m)
# des.posp.m<-svydesign(id=~1, strata=~ Strat_conc, weights=~sw,data=site.sw.posp.m)
# 
# des.d<-svydesign(id=~1, strata=~ Strat_conc, weights=~sw,data=site.sw.d)
# des.pocs.d<-svydesign(id=~1, strata=~ Strat_conc, weights=~sw,data=site.sw.pocs.d)
# des.mosp.d<-svydesign(id=~1, strata=~ Strat_conc, weights=~sw,data=site.sw.mosp.d)
# des.posp.d<-svydesign(id=~1, strata=~ Strat_conc, weights=~sw,data=site.sw.posp.d)


#  Test fixed effects of year
modR.s <- svyglm(AdColCount ~ OBS_YEAR, offset = TRANSECTAREA_ad, design = subset(des, DEPTH_BIN == "Shallow"), family = "poisson")
car::Anova(modR.s, type=3, test.statistic = "F") # NS
modR.m <- svyglm(AdColCount ~ OBS_YEAR, offset = TRANSECTAREA_ad, design = subset(des, DEPTH_BIN == "Mid"), family = "poisson")
car::Anova(modR.m, type=3, test.statistic = "F") # NS
modR.d <- svyglm(AdColCount ~ OBS_YEAR, offset = TRANSECTAREA_ad, design = subset(des, DEPTH_BIN == "Deep"), family = "poisson")
car::Anova(modR.d, type=3, test.statistic = "F") # Year = 0.0004

modR.s.pocs <- svyglm(AdColCount ~ OBS_YEAR, offset = TRANSECTAREA_ad, design = subset(des.pocs, DEPTH_BIN == "Shallow"), family = "poisson")
car::Anova(modR.s.pocs, type=3, test.statistic = "F") # Year = 0.0024
modR.m.pocs <- svyglm(AdColCount ~ OBS_YEAR, offset = TRANSECTAREA_ad, design = subset(des.pocs, DEPTH_BIN == "Mid"), family = "poisson")
car::Anova(modR.m.pocs, type=3, test.statistic = "F") # NS
modR.d.pocs <- svyglm(AdColCount ~ OBS_YEAR, offset = TRANSECTAREA_ad, design = subset(des.pocs, DEPTH_BIN == "Deep"), family = "poisson")
car::Anova(modR.d.pocs, type=3, test.statistic = "F") # Year = 0.0361

modR.s.mosp <- svyglm(AdColCount ~ OBS_YEAR, offset = TRANSECTAREA_ad, design = subset(des.mosp, DEPTH_BIN == "Shallow"), family = "poisson")
car::Anova(modR.s.mosp, type=3, test.statistic = "F") # NS
modR.m.mosp <- svyglm(AdColCount ~ OBS_YEAR, offset = TRANSECTAREA_ad, design = subset(des.mosp, DEPTH_BIN == "Mid"), family = "poisson")
car::Anova(modR.m.mosp, type=3, test.statistic = "F") # NS
modR.d.mosp <- svyglm(AdColCount ~ OBS_YEAR, offset = TRANSECTAREA_ad, design = subset(des.mosp, DEPTH_BIN == "Deep"), family = "poisson")
car::Anova(modR.d.mosp, type=3, test.statistic = "F") # Year = 0.0041

modR.s.posp <- svyglm(AdColCount ~ OBS_YEAR, offset = TRANSECTAREA_ad, design = subset(des.posp, DEPTH_BIN == "Shallow"), family = "poisson")
car::Anova(modR.s.posp, type=3, test.statistic = "F") # NS
modR.m.posp <- svyglm(AdColCount ~ OBS_YEAR, offset = TRANSECTAREA_ad, design = subset(des.posp, DEPTH_BIN == "Mid"), family = "poisson")
car::Anova(modR.m.posp, type=3, test.statistic = "F") # NS
modR.d.posp <- svyglm(AdColCount ~ OBS_YEAR, offset = TRANSECTAREA_ad, design = subset(des.posp, DEPTH_BIN == "Deep"), family = "poisson")
car::Anova(modR.d.posp, type=3, test.statistic = "F") # NS



with(subset(site.sw.pocs, DEPTH_BIN == "Deep"), tapply((AdColDen), OBS_YEAR, shapiro.test))
bartlett.test(AdColDen ~ Strat_conc, subset(site.sw.pocs, DEPTH_BIN == "Deep"))


# # Previous method of Svyglm
# modR<-svyglm(AdColCount ~ OBS_YEAR*DEPTH_BIN, offset = TRANSECTAREA_ad, family = "poisson", design=des) # all significant
# modR.pocs<-svyglm(AdColCount ~ OBS_YEAR*DEPTH_BIN, offset = TRANSECTAREA_ad,family = "poisson", design=des.pocs) # Year and depth significant
# modR.mosp<-svyglm(AdColCount ~ OBS_YEAR*DEPTH_BIN, offset = TRANSECTAREA_ad,family = "poisson", design=des.mosp) # Depth and year:depth significant
# modR.posp<-svyglm(AdColCount ~ OBS_YEAR*DEPTH_BIN, offset = TRANSECTAREA_ad,family = "poisson", design=des.posp) # Depth significant
# 
# tot <- anova(modR)
# poc <- anova(modR.pocs)
# mos <- anova(modR.mosp)
# pos <- anova(modR.posp)
# summary(modR); summary(modR.pocs); summary(modR.mosp); summary(modR.posp)



# Post-hoc
tot.d <- emmeans(modR.d, specs = pairwise~OBS_YEAR, adjust = "none") # 2015 is sig less than 2018 and 2023
tot.d.df <- as.data.frame(tot.d[[2]])

pocs.s <- emmeans(modR.s.pocs, specs = pairwise~OBS_YEAR, adjust = "none") # all years sig different. 2015>2023>2018
pocs.s.df <- as.data.frame(pocs.s[[2]])
pocs.d <- emmeans(modR.d.pocs, specs = pairwise~OBS_YEAR, adjust = "none") # 2015 is sig less than 2023
pocs.d.df <- as.data.frame(pocs.d[[2]])

mosp.d <- emmeans(modR.d.mosp, specs = pairwise~OBS_YEAR, adjust = "none") # 2015 is sig less than 2018 and 2023 
mosp.d.df <- as.data.frame(mosp.d[[2]])

# Previous version: Year
# y <- svyglm(AdColCount ~ OBS_YEAR, offset = TRANSECTAREA_ad, family = "poisson", design = des) # All corals
# summary(glht(y, mcp(OBS_YEAR="Tukey")))  # 2018 sig different than 2015 and 2023
# 
# 
# y <- svyglm(AdColCount ~ OBS_YEAR, offset = TRANSECTAREA_ad, family = "poisson", design = des.pocs) # POCS
# summary(glht(y, mcp(OBS_YEAR="Tukey")))  # 2018 sig different than 2015 and 2023
# 
# y <- svyglm(AdColCount ~ OBS_YEAR, offset = TRANSECTAREA_ad, family = "poisson", design = des.posp) # POSP
# summary(glht(y, mcp(OBS_YEAR="Tukey"))) # None significant
#
# 
# y <- svyglm(AdColCount ~ OBS_YEAR, offset = TRANSECTAREA_ad, family = "poisson", design = des.mosp) # MOSP
# summary(glht(y, mcp(OBS_YEAR="Tukey"))) # None significant
# 
# summary(glht(modR.pocs, mcp(OBS_YEAR="Tukey"))) # 2023 sig different from 2015
# summary(glht(modR.mosp, mcp(OBS_YEAR="Tukey"))) # 2023 and 2018 sig different than 2015
# summary(glht(modR.posp, mcp(OBS_YEAR="Tukey"))) # 2023 and 2018 sig different than 2015


# Apply test corrections
tot.d.df$p.adj <- p.adjust(tot.d.df$p.value, method = "BH") # above holds true
pocs.s.df$p.adj <- p.adjust(pocs.s.df$p.value, method = "BH"); pocs.s.df # above holds true
pocs.d.df$p.adj <- p.adjust(pocs.d.df$p.value, method = "BH"); pocs.d.df # above holds true
mosp.d.df$p.adj <- p.adjust(mosp.d.df$p.value, method = "BH"); mosp.d.df # above holds true



#Calculate regional mean and SE
# depth_mean<-svyby(~CORAL,~OBS_YEAR+DEPTH_BIN,des,svymean) # Dont have "CORAL" column
depth_mean<-svyby(~AdColDen,~OBS_YEAR+DEPTH_BIN,des,svymean)
depth_mean_pocs<-svyby(~AdColDen,~OBS_YEAR+DEPTH_BIN,des.pocs,svymean)
depth_mean_mosp<-svyby(~AdColDen,~OBS_YEAR+DEPTH_BIN,des.mosp,svymean)
depth_mean_posp<-svyby(~AdColDen,~OBS_YEAR+DEPTH_BIN,des.posp,svymean)
# depth_mean<-svyby(~AdColCount,~OBS_YEAR+DEPTH_BIN,des,svymean)
depth_mean


# Put all depth_mean dataframes together for plotting
depth_mean$Genus <- paste("Total")
depth_mean$sig <- c("a","b","b",NA,NA,NA,NA,NA,NA)
depth_mean_pocs$Genus <- paste("Pocillopora")
depth_mean_pocs$sig <- c("a","ab","b",NA,NA,NA,"a","b","c")
depth_mean_mosp$Genus <- paste("Montipora")
depth_mean_mosp$sig <- c("a","b","b",NA,NA,NA,NA,NA,NA)
depth_mean_posp$Genus <- paste("Porites")
depth_mean_posp$sig <- c(NA,NA,NA,NA,NA,NA,NA,NA,NA)
depth_mean_tot <- rbind(depth_mean, depth_mean_pocs, depth_mean_mosp, depth_mean_posp)

# Plot Data
s <- ggplot(depth_mean_tot %>% filter(DEPTH_BIN == "Shallow"), 
            aes(x = OBS_YEAR, y = AdColDen, group = Genus, fill = OBS_YEAR)) +
        geom_bar(stat = "identity") + 
        scale_fill_manual(values = alpha(c("#440154FF","#22A884FF","#FDE725FF"))) +
        geom_errorbar(data = depth_mean_tot %>% filter(DEPTH_BIN == "Shallow"), 
                      aes(ymin = AdColDen-se, ymax = AdColDen+se),
                      width = .2) +
        geom_text(aes(x=OBS_YEAR,y=AdColDen+se,label=sig, group = Genus),
                  position = position_dodge(),
                  vjust = -0.5) +
        facet_wrap(~Genus, nrow = 1) +
        guides(fill="none") +
        theme_bw() +
        theme(panel.grid.minor = element_blank(),
              panel.grid.major = element_blank(),
              axis.title.y = element_blank(),
              axis.title.x = element_blank()) +
        ggtitle("Shallow")

m <- ggplot(depth_mean_tot %>% filter(DEPTH_BIN == "Mid"), 
            aes(x = OBS_YEAR, y = AdColDen, group = Genus, fill = OBS_YEAR)) +
        geom_bar(stat = "identity") + 
        scale_fill_manual(values = alpha(c("#440154FF","#22A884FF","#FDE725FF"))) +
        geom_errorbar(data = depth_mean_tot %>% filter(DEPTH_BIN == "Mid"),
                      aes(ymin = AdColDen-se, ymax = AdColDen+se),
                      width = .1) +
        geom_text(aes(x=OBS_YEAR,y=AdColDen+se,label=sig, group = Genus),
                  position = position_dodge(),
                  vjust = -0.5) +
        facet_wrap(~Genus, nrow = 1) +
        guides(fill="none") +
        theme_bw() +
        theme(panel.grid.minor = element_blank(),
              panel.grid.major = element_blank(),
              axis.title.x = element_blank(),
              axis.title.y = element_text(size=12, face = "bold"))  +
        ggtitle("Mid")

d <- ggplot(depth_mean_tot %>% filter(DEPTH_BIN == "Deep"), 
            aes(x = OBS_YEAR, y = AdColDen, group = Genus, fill = OBS_YEAR)) +
        geom_bar(stat = "identity") + 
        scale_fill_manual(values = alpha(c("#440154FF","#22A884FF","#FDE725FF"))) +
        geom_errorbar(data = depth_mean_tot %>% filter(DEPTH_BIN == "Deep"),
                      aes(ymin = AdColDen-se, ymax = AdColDen+se),
                      width = .1) +
        geom_text(aes(x=OBS_YEAR,y=AdColDen+se,label=sig, group = Genus),
                 position = position_dodge(),
                 vjust = -0.5) +
        facet_wrap(~Genus, nrow = 1) +
        guides(fill="none") +
        ylim(0,22) +
        theme_bw() +
        theme(panel.grid.minor = element_blank(),
              panel.grid.major = element_blank(),
              axis.title.y = element_blank(),
              axis.title.x = element_text(size = 12, face = "bold")) +
        ggtitle("Deep")

# tot <- ggplot(depth_mean, 
#               aes(x = OBS_YEAR, y = AdColDen, fill = OBS_YEAR)) +
#   geom_bar(stat = "identity") + 
#   scale_fill_manual(values = alpha(c("#440154FF","#22A884FF","#FDE725FF"))) +
#   geom_errorbar(data = depth_mean,
#                 aes(ymin = AdColDen-se, ymax = AdColDen+se),
#                 width = .1) +
#   facet_wrap(~DEPTH_BIN, nrow = 1) +
#   guides(fill="none") +
#   theme_bw() +
#   theme(panel.grid.minor = element_blank(),
#         panel.grid.major = element_blank()) 


ggarrange(s, m, d, nrow = 3, align = "v", heights = c(3.6, 3.6,4))
         


#### Genus level Partial Mortality####


#Subset data
site.sw$OBS_YEAR<-as.factor(site.sw$OBS_YEAR)
site.sw$DEPTH_BIN<-as.factor(site.sw$DEPTH_BIN)
site.sw[is.na(site.sw)] <- 0
site.sw$dead <- site.sw$Ave.od + site.sw$Ave.rd

site.sw.pocs <- site.sw %>% filter(GENUS_CODE == "POCS") # above code isn't 50 so here I'm looking at just total coral
site.sw.mosp <- site.sw %>% filter(GENUS_CODE == "MOSP") # 
site.sw.posp <- site.sw %>% filter(GENUS_CODE == "POSP") # 
site.sw <- site.sw %>% filter(GENUS_CODE == "SSSS")


#Establish survey design with survey weights
des<-svydesign(id=~1, strata=~ Strat_conc, weights=~sw,data=site.sw)
des.pocs<-svydesign(id=~1, strata=~ Strat_conc, weights=~sw,data=site.sw.pocs)
des.mosp<-svydesign(id=~1, strata=~ Strat_conc, weights=~sw,data=site.sw.mosp)
des.posp<-svydesign(id=~1, strata=~ Strat_conc, weights=~sw,data=site.sw.posp)



#  Test fixed effects of year
modR.s <- svyglm(dead ~ OBS_YEAR, design = subset(des, DEPTH_BIN == "Shallow"), family = "poisson")
car::Anova(modR.s, type=3, test.statistic = "F") # NS
modR.m <- svyglm(dead ~ OBS_YEAR,  design = subset(des, DEPTH_BIN == "Mid"), family = "poisson")
car::Anova(modR.m, type=3, test.statistic = "F") # NS
modR.d <- svyglm(dead ~ OBS_YEAR,  design = subset(des, DEPTH_BIN == "Deep"), family = "poisson")
car::Anova(modR.d, type=3, test.statistic = "F") # NS

modR.s.pocs <- svyglm(dead ~ OBS_YEAR, design = subset(des.pocs, DEPTH_BIN == "Shallow"), family = "poisson")
car::Anova(modR.s.pocs, type=3, test.statistic = "F") # Year <0.001
modR.m.pocs <- svyglm(dead ~ OBS_YEAR,  design = subset(des.pocs, DEPTH_BIN == "Mid"), family = "poisson")
car::Anova(modR.m.pocs, type=3, test.statistic = "F") # NS
modR.d.pocs <- svyglm(dead ~ OBS_YEAR, design = subset(des.pocs, DEPTH_BIN == "Deep"), family = "poisson")
car::Anova(modR.d.pocs, type=3, test.statistic = "F") # NS

modR.s.mosp <- svyglm(dead ~ OBS_YEAR,design = subset(des.mosp, DEPTH_BIN == "Shallow"), family = "poisson")
car::Anova(modR.s.mosp, type=3, test.statistic = "F") # NS
modR.m.mosp <- svyglm(dead ~ OBS_YEAR,design = subset(des.mosp, DEPTH_BIN == "Mid"), family = "poisson")
car::Anova(modR.m.mosp, type=3, test.statistic = "F") # NS
modR.d.mosp <- svyglm(dead ~ OBS_YEAR,  design = subset(des.mosp, DEPTH_BIN == "Deep"), family = "poisson")
car::Anova(modR.d.mosp, type=3, test.statistic = "F") # NS

modR.s.posp <- svyglm(dead ~ OBS_YEAR, design = subset(des.posp, DEPTH_BIN == "Shallow"), family = "poisson")
car::Anova(modR.s.posp, type=3, test.statistic = "F") # NS
modR.m.posp <- svyglm(dead ~ OBS_YEAR,  design = subset(des.posp, DEPTH_BIN == "Mid"), family = "poisson")
car::Anova(modR.m.posp, type=3, test.statistic = "F") # NS
modR.d.posp <- svyglm(dead ~ OBS_YEAR,  design = subset(des.posp, DEPTH_BIN == "Deep"), family = "poisson")
car::Anova(modR.d.posp, type=3, test.statistic = "F") # NS



# Post-hoc
pocs.s <- emmeans(modR.s.pocs, specs = pairwise~OBS_YEAR, adjust = "none"); pocs.s # 2018 is sig more than 2015 and 2023
pocs.s.df <- as.data.frame(pocs.s[[2]])

pocs.d <- emmeans(modR.d.pocs, specs = pairwise~OBS_YEAR, adjust = "none"); pocs.d # 2018 > 2015
pocs.d.df <- as.data.frame(pocs.d[[2]])

# Apply test corrections
pocs.s.df$p.adj <- p.adjust(pocs.s.df$p.value, method = "BH"); pocs.s.df # above holds true
pocs.d.df$p.adj <- p.adjust(pocs.d.df$p.value, method = "BH"); pocs.d.df # above holds true


#Calculate regional mean and SE
od_mean<-svyby(~dead,~OBS_YEAR+DEPTH_BIN,des,svymean)
od_mean_pocs<-svyby(~dead,~OBS_YEAR+DEPTH_BIN,des.pocs,svymean)
od_mean_mosp<-svyby(~dead,~OBS_YEAR+DEPTH_BIN,des.mosp,svymean)
od_mean_posp<-svyby(~dead,~OBS_YEAR+DEPTH_BIN,des.posp,svymean)
od_mean


# Put all od_mean dataframes together for plotting
od_mean$Genus <- paste("Total")
od_mean$sig <- c(NA,NA,NA,NA,NA,NA,NA,NA,NA)
od_mean_pocs$Genus <- paste("Pocillopora")
od_mean_pocs$sig <- c(NA,NA,NA,NA,NA,NA,"a","b","a")
od_mean_mosp$Genus <- paste("Montipora")
od_mean_mosp$sig <- c(NA,NA,NA,NA,NA,NA,NA,NA,NA)
od_mean_posp$Genus <- paste("Porites")
od_mean_posp$sig <- c(NA,NA,NA,NA,NA,NA,NA,NA,NA)
od_mean_tot <- rbind(od_mean, od_mean_pocs, od_mean_mosp, od_mean_posp)


# od_mean_tot$DEPTH_BIN <- factor(od_mean_tot$DEPTH_BIN, levels = c("Shallow","Mid","Deep")) # Add Posthoc groupings from glms
# od_mean_tot <- od_mean[order(od_mean_tot$DEPTH_BIN),];od_mean_tot 


# Plot Partial mortality data

hist(ssss$dead)

s <- ggplot(od_mean_tot %>% filter(DEPTH_BIN == "Shallow"), 
            aes(x=OBS_YEAR, y=dead, fill=OBS_YEAR, group = Genus)) +
            geom_bar(stat = "identity") +
            scale_fill_manual(values = alpha(c("#440154FF","#22A884FF","#FDE725FF"))) +
            geom_errorbar(data = od_mean_tot %>% filter(DEPTH_BIN == "Shallow"),
                          aes(ymin=dead-se, ymax=dead+se), width = .1) +
            geom_text(aes(x=OBS_YEAR,y=dead+se+3,label=sig, group = Genus),
                      position = position_dodge()) +
            facet_wrap(~Genus, nrow = 1) +
            theme_bw() +
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.spacing = unit(0, "lines"),
                  strip.background = element_blank(),
                  strip.placement = "outside",
                  strip.text = element_text(size = 12),
                  legend.position = "none",
                  axis.line = element_line(color = "black"),
                  text = element_text(size = 12),
                  axis.text.y = element_text(colour="black"),
                  axis.text.x = element_text(size=12),
                  axis.title.y = element_text(size = 14, face = "bold")) +
            scale_fill_viridis_d() +
            labs(x="",y="")+
            scale_y_continuous(expand = c(0,0), limits = c(0,41))+
            ggtitle("Shallow")

m <- ggplot(od_mean_tot %>% filter(DEPTH_BIN == "Mid"), 
            aes(x=OBS_YEAR, y=dead, fill=OBS_YEAR, group = Genus)) +
      geom_bar(stat = "identity") +
      scale_fill_manual(values = alpha(c("#440154FF","#22A884FF","#FDE725FF"))) +
      geom_errorbar(data = od_mean_tot %>% filter(DEPTH_BIN == "Mid"),
                    aes(ymin=dead-se, ymax=dead+se), width = .1) +
      geom_text(aes(x=OBS_YEAR,y=dead+se+.8,label=sig, group = Genus),
                position = position_dodge()) +
      facet_wrap(~Genus, nrow = 1) +
      theme_bw() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.spacing = unit(0, "lines"),
            strip.background = element_blank(),
            strip.placement = "outside",
            strip.text = element_text(size = 12),
            legend.position = "none",
            axis.line = element_line(color = "black"),
            text = element_text(size = 12),
            axis.text.y = element_text(colour="black"),
            axis.text.x = element_text(size=12),
            axis.title.y = element_text(size = 13, face = "bold")) +
      scale_fill_viridis_d() +
      labs(x="",y=expression(paste("Mean % Old Dead")))+
      scale_y_continuous(expand = c(0,0), limits = c(0,41))+
      ggtitle("Mid")

d <- ggplot(od_mean_tot %>% filter(DEPTH_BIN == "Deep"), 
            aes(x=OBS_YEAR, y=dead, fill=OBS_YEAR, group = Genus)) +
      geom_bar(stat = "identity") +
      scale_fill_manual(values = alpha(c("#440154FF","#22A884FF","#FDE725FF"))) +
      geom_errorbar(data = od_mean_tot %>% filter(DEPTH_BIN == "Deep"),
                    aes(ymin=dead-se, ymax=dead+se), width = .1) +
      geom_text(aes(x=OBS_YEAR,y=dead+se+.8,label=sig, group = Genus),
                position = position_dodge()) +
      facet_wrap(~Genus, nrow = 1) +
      theme_bw() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.spacing = unit(0, "lines"),
            strip.background = element_blank(),
            strip.placement = "outside",
            strip.text = element_text(size = 12),
            legend.position = "none",
            axis.line = element_line(color = "black"),
            text = element_text(size = 12),
            axis.text.y = element_text(colour="black"),
            axis.text.x = element_text(size=12),
            axis.title.x = element_text(size = 13, face = "bold")) +
      scale_fill_viridis_d() +
      labs(x="Year",y="")+
      scale_y_continuous(expand = c(0,0), limits = c(0,41))+
      ggtitle("Deep")


ggarrange(s, m, d, nrow = 3, align = "v") 
          heights = c(3.6, 3.6,4)
          
          
          
#### Adult density using SIMPER results ####

# Merge demo data with new NH values pooled across the 2 swains sectors
swan<-left_join(swan,NH) 
swan <- filter(swan, DEPTH_BIN == "Shallow")

# Calculate survey weights (inverse proportion weighting)
w.df<-swan %>%
  group_by(OBS_YEAR,DEPTH_BIN,NH) %>%
  summarise(n = length(unique(SITE)))

w.df$sw<-w.df$NH/w.df$n #calculate survey weights for each site

site.sw<-left_join(swan,w.df) #merge weights with site-level data
head(site.sw)

#Create contactenated Strata variable 
site.sw$Strat_conc<-paste(site.sw$OBS_YEAR, site.sw$DEPTH_BIN,sep = "_")
site.sw$OBS_YEAR<-as.factor(site.sw$OBS_YEAR)

#Subset data
site.sw.mospfo <- site.sw %>% filter(TAXONCODE_2 == "MOSP_FO") 
site.sw.mospem <- site.sw %>% filter(TAXONCODE_2 == "MOSP_EM") # 
site.sw.pmvc <- site.sw %>% filter(TAXONCODE_2 == "PMVC") # 
site.sw.pgwc <- site.sw %>% filter(TAXONCODE_2 == "PGWC") # 
site.sw.pospem <- site.sw %>% filter(TAXONCODE_2 == "POSP_EM") # 
site.sw.pospmd <- site.sw %>% filter(TAXONCODE_2 == "POSP_MD")

site.sw.pgwc$OBS_YEAR<-as.factor(site.sw.pmvc$OBS_YEAR)
site.sw.pmvc$OBS_YEAR<-as.factor(site.sw.pgwc$OBS_YEAR)
site.sw.mospfo$OBS_YEAR<-as.factor(site.sw.mospfo$OBS_YEAR)
site.sw.mospem$OBS_YEAR<-as.factor(site.sw.mospem$OBS_YEAR)
site.sw.pospem$OBS_YEAR<-as.factor(site.sw.pospem$OBS_YEAR)
site.sw.pospmd$OBS_YEAR<-as.factor(site.sw.pospmd$OBS_YEAR)


#remove strata that have less than 2 sites
site.sw.pgwc<-subset(site.sw.pgwc,n>1)
summary(site.sw.pgwc$n)
length(unique(site.sw$SITE))

site.sw.pgwc<-subset(site.sw.pgwc,n>1)

#Check for duplicate sites- df should be empty
site.sw.pgwc %>% 
  group_by(SITE) %>% 
  filter(n()>1)


# Test assumptions for parametric tests
with(site.sw.pmvc, tapply((AdColDen), OBS_YEAR, shapiro.test)) # NS
with(site.sw.pgwc, tapply((AdColDen), OBS_YEAR, shapiro.test)) # not working (0 in 2018)
with(site.sw.mospfo, tapply((AdColDen), OBS_YEAR, shapiro.test)) # NS
with(site.sw.mospem, tapply((AdColDen), OBS_YEAR, shapiro.test)) # NS
with(site.sw.pospem, tapply((AdColDen), OBS_YEAR, shapiro.test)) # not working
with(site.sw.pospmd, tapply((AdColDen), OBS_YEAR, shapiro.test)) # not working (0 in 2018)
bartlett.test(AdColDen ~ Strat_conc,site.sw.pgwc)


#Establish survey design with survey weights
des.pgwc<-svydesign(id=~1, strata=~ Strat_conc, weights=~sw,data=site.sw.pgwc)
des.pmvc<-svydesign(id=~1, strata=~ Strat_conc, weights=~sw,data=site.sw.pmvc)
des.mospfo<-svydesign(id=~1, strata=~ Strat_conc, weights=~sw,data=site.sw.mospfo)
des.mospem<-svydesign(id=~1, strata=~ Strat_conc, weights=~sw,data=site.sw.mospem)
des.pospem<-svydesign(id=~1, strata=~ Strat_conc, weights=~sw,data=site.sw.pospem)
des.pospmd<-svydesign(id=~1, strata=~ Strat_conc, weights=~sw,data=site.sw.pospmd)


# Non-parametric version of svyglm
svyranktest(AdColDen ~ OBS_YEAR, design=subset(des.pgwc,DEPTH_BIN=="Shallow"), test=("KruskalWallis")) # p < 0.002
# svyranktest(AdColDen ~ OBS_YEAR, design=subset(des.pmvc,DEPTH_BIN=="Shallow"), test=("KruskalWallis")) # p < 0.001
# svyranktest(AdColDen ~ OBS_YEAR, design=subset(des.mospem,DEPTH_BIN=="Shallow"), test=("KruskalWallis")) # NS
# svyranktest(AdColDen ~ OBS_YEAR, design=subset(des.mospfo,DEPTH_BIN=="Shallow"), test=("KruskalWallis")) # NS
svyranktest(AdColDen ~ OBS_YEAR, design=subset(des.pospem,DEPTH_BIN=="Shallow"), test=("KruskalWallis")) # NS
svyranktest(AdColDen ~ OBS_YEAR, design=subset(des.pospmd,DEPTH_BIN=="Shallow"), test=("KruskalWallis")) # p < 0.001


#  Test fixed effects of year
modR.pgwc <- svyglm(AdColCount ~ OBS_YEAR, offset = TRANSECTAREA_ad, design = des.pgwc, family = "poisson")

modR.pmvc <- svyglm(AdColCount ~ OBS_YEAR,  offset = TRANSECTAREA_ad, design = des.pmvc, family = "poisson")
car::Anova(modR.pmvc, type=3, test.statistic = "F") # p = 0.03
modR.pmvc <- svyglm(AdColDen ~ OBS_YEAR,  design = des.pmvc)
car::Anova(modR.pmvc, type=3, test.statistic = "F") # p = 0.03

modR.mospem <- svyglm(AdColCount ~ OBS_YEAR, offset = TRANSECTAREA_ad,   design = des.mospem, family = "poisson")
car::Anova(modR.mospem, type=3, test.statistic = "F") # NS

modR.mospfo <- svyglm(AdColCount ~ OBS_YEAR, offset = TRANSECTAREA_ad, design = des.mospfo, family = "poisson")
car::Anova(modR.mospfo, type=3, test.statistic = "F") # NS

# modR.pospem <- svyglm(AdColCount ~ OBS_YEAR, offset = TRANSECTAREA_ad, design = des.pospem, family = "poisson")
# car::Anova(modR.pospem, type=3, test.statistic = "F") # NS
# 
# modR.pospmd <- svyglm(AdColCount ~ OBS_YEAR, offset = TRANSECTAREA_ad, design = des.pospmd, family = "poisson")
# car::Anova(modR.pospmd, type=3, test.statistic = "F") # NS


# Post-hoc
tot.pmvc <- emmeans(modR.pmvc, specs = pairwise~OBS_YEAR, adjust = "none") # 2015 is sig less than 2018 and 2023
tot.pmvc <- emmeans(modR.pmvc,pairwise~OBS_YEAR) 
tot.pmvc.df <- as.data.frame(tot.pmvc[[2]])



# Apply test corrections
tot.pmvc.df$p.adj <- p.adjust(tot.pmvc.df$p.value, method = "BH"); tot.pmvc.df # 2018 no longer smaller than 2023


l<-c("2015","2018","2023") # p value for non-parametric tests
ps<-matrix(NA,3,3)
dimnames(ps)<-list(l,l)
for(i in 1:2){
  for(j in (i+1):3){
    ps[i,j]<-svyranktest(AdColDen~OBS_YEAR, subset(des.pospmd, OBS_YEAR %in% l[c(i,j)]))$p.value
  }
}
ps

p.adjust(ps[!is.na(ps)],method="hochberg")

# PGWC: all years sig different
# POSPMD: 2018 sig different


#Calculate regional mean and SE
depth_mean_pmvc<-svyby(~AdColDen,~OBS_YEAR+DEPTH_BIN,des.pmvc,svymean)
depth_mean_pgwc<-svyby(~AdColDen,~OBS_YEAR+DEPTH_BIN,des.pgwc,svymean)
depth_mean_pospmd<-svyby(~AdColDen,~OBS_YEAR+DEPTH_BIN,des.pospmd,svymean)
depth_mean_pospem<-svyby(~AdColDen,~OBS_YEAR+DEPTH_BIN,des.pospem,svymean) 
depth_mean_mospfo<-svyby(~AdColDen,~OBS_YEAR+DEPTH_BIN,des.mospfo,svymean)
depth_mean_mospem<-svyby(~AdColDen,~OBS_YEAR+DEPTH_BIN,des.mospem,svymean) # is parametric and NS, but becomes significant if no parametric glm is run


# Put all depth_mean dataframes together for plotting
depth_mean_pmvc$Genus <- paste("P. meandrina/verrucosa complex")
depth_mean_pmvc$sig <- c("a","b","ab")
depth_mean_pgwc$Genus <- paste("P. grandis/woodjonsei complex")
depth_mean_pgwc$sig <- c("a","b","c")
depth_mean_pospmd$Genus <- paste("Mounding Porites spp.")
depth_mean_pospmd$sig <- c("a","b","a")
depth_mean_pospem$Genus <- paste("Encrusting Porites spp.")
depth_mean_pospem$sig <- c(NA,NA,NA)
depth_mean_mospem$Genus <- paste("Encrusting Montipora spp.")
depth_mean_mospem$sig <- c(NA,NA,NA)
depth_mean_mospfo$Genus <- paste("Foliose Montipora spp.")
depth_mean_mospfo$sig <- c(NA,NA,NA)
depth_mean_tot <- rbind(depth_mean_pgwc,depth_mean_pospem,depth_mean_pospmd,
                        depth_mean_mospem,depth_mean_mospfo,depth_mean_pmvc)

# Plot Data
s <- ggplot(depth_mean_tot, aes(x = OBS_YEAR, y = AdColDen, group = Genus, fill = OBS_YEAR)) +
        geom_bar(stat = "identity") + 
        scale_fill_manual(values = alpha(c("#440154FF","#22A884FF","#FDE725FF"))) +
        geom_errorbar(data = depth_mean_tot %>% filter(DEPTH_BIN == "Shallow"), 
                      aes(ymin = AdColDen-se, ymax = AdColDen+se),
                      width = .3) +
        geom_text(aes(x=OBS_YEAR,y=AdColDen+se,label=sig, group = Genus),
                  position = position_dodge(), size = 6,
                  vjust = -0.5) +
        facet_wrap(~Genus,nrow = 2, labeller = label_wrap_gen(width = 16)) +
        guides(fill="none") +
        theme_bw() +
        ylim(0,8.9) +
        theme(panel.grid.minor = element_blank(),
              panel.grid.major = element_blank(),
              strip.text.x = element_text(size = 15, face = "bold"),
              axis.text.x = element_text(size = 13, face = "bold"),
              axis.text.y = element_text(size = 13, face = "bold"),
              axis.title.x = element_text(size=15, face = "bold"),
              axis.title.y = element_text(size=15, face = "bold"))+
        ylab("Mean Adult Colonies m - 2") +
        xlab("Year")
        # ggtitle("Shallow")
        
        # levels = c("P. meandrina/verrucosa complex","Mounding Porites spp.", "Foliose Montipora spp.",
        #            "P. grandis/woodjonsei complex","Encrusting Porites spp.", "Encrusting Montipora spp.")),




#### Partial Mortality using SIMPER results ####

# Merge demo data with new NH values pooled across the 2 swains sectors
swam<-left_join(swam,NH) 
swam <- filter(swam, DEPTH_BIN == "Shallow")

# Calculate survey weights (inverse proportion weighting)
w.df<-swam %>%
  group_by(OBS_YEAR,DEPTH_BIN,NH) %>%
  summarise(n = length(unique(SITE)))

w.df$sw<-w.df$NH/w.df$n #calculate survey weights for each site

site.sw<-left_join(swam,w.df) #merge weights with site-level data
head(site.sw)

#Create contactenated Strata variable 
site.sw$Strat_conc<-paste(site.sw$OBS_YEAR, site.sw$DEPTH_BIN,sep = "_")
site.sw$OBS_YEAR<-as.factor(site.sw$OBS_YEAR)

#Subset data
site.sw.mospfo <- site.sw %>% filter(TAXONCODE_2 == "MOSP_FO") 
site.sw.mospem <- site.sw %>% filter(TAXONCODE_2 == "MOSP_EM") # 
site.sw.pmvc <- site.sw %>% filter(TAXONCODE_2 == "PMVC") # 
site.sw.pgwc <- site.sw %>% filter(TAXONCODE_2 == "PGWC") # 
site.sw.pospem <- site.sw %>% filter(TAXONCODE_2 == "POSP_EM") # 
site.sw.pospmd <- site.sw %>% filter(TAXONCODE_2 == "POSP_MD")


#Establish survey design with survey weights
des.pgwc<-svydesign(id=~1, strata=~ Strat_conc, weights=~sw,data=site.sw.pgwc)
des.pmvc<-svydesign(id=~1, strata=~ Strat_conc, weights=~sw,data=site.sw.pmvc)
des.mospfo<-svydesign(id=~1, strata=~ Strat_conc, weights=~sw,data=site.sw.mospfo)
des.mospem<-svydesign(id=~1, strata=~ Strat_conc, weights=~sw,data=site.sw.mospem)
des.pospem<-svydesign(id=~1, strata=~ Strat_conc, weights=~sw,data=site.sw.pospem)
des.pospmd<-svydesign(id=~1, strata=~ Strat_conc, weights=~sw,data=site.sw.pospmd)


# Test assumptions
with(site.sw.pmvc, tapply((Ave.od), OBS_YEAR, shapiro.test)) # passed
bartlett.test(Ave.od ~ Strat_conc,site.sw.pmvc) 
with(site.sw.pgwc, tapply((Ave.od), OBS_YEAR, shapiro.test)) # not passed
bartlett.test(Ave.od ~ Strat_conc,site.sw.pgwc) 
with(site.sw.mospem, tapply((Ave.od), OBS_YEAR, shapiro.test)) # passed
bartlett.test(Ave.od ~ Strat_conc,site.sw.mospem)
with(site.sw.mospfo, tapply((Ave.od), OBS_YEAR, shapiro.test)) # passed
bartlett.test(Ave.od ~ Strat_conc,site.sw.mospfo)
with(site.sw.pospem, tapply((Ave.od), OBS_YEAR, shapiro.test)) # not passed
bartlett.test(Ave.od ~ Strat_conc,site.sw.pospem)
with(site.sw.pospmd, tapply((Ave.od), OBS_YEAR, shapiro.test)) # not passed
bartlett.test(Ave.od ~ Strat_conc,site.sw.pospmd)


# Non-parametric version of svyglm
svyranktest(Ave.od ~ OBS_YEAR, design=subset(des.pgwc,DEPTH_BIN=="Shallow"), test=("KruskalWallis")) # NS
# svyranktest(Ave.od ~ OBS_YEAR, design=subset(des.pmvc,DEPTH_BIN=="Shallow"), test=("KruskalWallis")) # p < 0.001
# svyranktest(Ave.od ~ OBS_YEAR, design=subset(des.mospem,DEPTH_BIN=="Shallow"), test=("KruskalWallis")) # NS
# svyranktest(Ave.od ~ OBS_YEAR, design=subset(des.mospfo,DEPTH_BIN=="Shallow"), test=("KruskalWallis")) # NS
svyranktest(Ave.od ~ OBS_YEAR, design=subset(des.pospem,DEPTH_BIN=="Shallow"), test=("KruskalWallis")) # NS
svyranktest(Ave.od ~ OBS_YEAR, design=des.pospmd, test=("KruskalWallis")) # p = 0.052


#  Test fixed effects of year
# modR.pgwc <- svyglm(Ave.od ~ OBS_YEAR,  design = des.pgwc)
# car::Anova(modR.pgwc, type=3, test.statistic = "F") # NS

modR.pmvc <- svyglm(Ave.od ~ OBS_YEAR,  design = des.pmvc)
car::Anova(modR.pmvc, type=3, test.statistic = "F") # < 0.001

modR.mospem <- svyglm(Ave.od ~ OBS_YEAR,  design = des.mospem)
car::Anova(modR.mospem, type=3, test.statistic = "F") # NS

modR.mospfo <- svyglm(Ave.od ~ OBS_YEAR, design = des.mospfo)
car::Anova(modR.mospfo, type=3, test.statistic = "F") # NS

# modR.pospem <- svyglm(Ave.od ~ OBS_YEAR, design = des.pospem)
# car::Anova(modR.pospem, type=3, test.statistic = "F") # NS
# 
# modR.pospmd <- svyglm(Ave.od ~ OBS_YEAR, design = des.pospmd)
# car::Anova(modR.pospmd, type=3, test.statistic = "F") # NS


# Post-hoc
# tot.pgwc <- emmeans(modR.pgwc, specs = pairwise~OBS_YEAR, adjust = "none") # 2015 is sig less than 2018 and 2023
# tot.pgwc.df <- as.data.frame(tot.pgwc[[2]])

tot.pmvc <- emmeans(modR.pmvc, specs = pairwise~OBS_YEAR, adjust = "none") # 2015 is sig less than 2018 and 2023
tot.pmvc.df <- as.data.frame(tot.pmvc[[2]])

# tot.pospmd <- emmeans(modR.pospmd, specs = pairwise~OBS_YEAR, adjust = "none") # 2015 is sig less than 2018 and 2023
# tot.pospmd.df <- as.data.frame(tot.pospmd[[2]])
# 
# tot.mospem <- emmeans(modR.mospem, specs = pairwise~OBS_YEAR, adjust = "none") # 2015 is sig less than 2018 and 2023
# tot.mospem.df <- as.data.frame(tot.mospem[[2]])


# Apply test corrections
tot.pmvc.df$p.adj <- p.adjust(tot.pmvc.df$p.value, method = "BH"); tot.pmvc.df # above holds true


#Calculate regional mean and SE
od_mean_pmvc<-svyby(~Ave.od,~OBS_YEAR+DEPTH_BIN,des.pmvc,svymean)
od_mean_pgwc<-svyby(~Ave.od,~OBS_YEAR+DEPTH_BIN,des.pgwc,svymean)
od_mean_pospmd<-svyby(~Ave.od,~OBS_YEAR+DEPTH_BIN,des.pospmd,svymean)
od_mean_pospem<-svyby(~Ave.od,~OBS_YEAR+DEPTH_BIN,des.pospem,svymean) 
od_mean_mospfo<-svyby(~Ave.od,~OBS_YEAR+DEPTH_BIN,des.mospfo,svymean)
od_mean_mospem<-svyby(~Ave.od,~OBS_YEAR+DEPTH_BIN,des.mospem,svymean) 

# Put all od_mean dataframes together for plotting
od_mean_pmvc$Genus <- paste("P. meandrina/verrucosa complex")
od_mean_pmvc$sig <- c("a","b","a")
od_mean_pgwc$Genus <- paste("P. grandis/woodjonsei complex")
od_mean_pgwc$sig <- c(NA,NA,NA)
od_mean_pospmd$Genus <- paste("Mounding Porites spp.")
od_mean_pospmd$sig <- c(NA,NA,NA)
od_mean_pospem$Genus <- paste("Encrusting Porites spp.")
od_mean_pospem$sig <- c(NA,NA,NA)
od_mean_mospem$Genus <- paste("Encrusting Montipora spp.")
od_mean_mospem$sig <- c(NA,NA,NA)
od_mean_mospfo$Genus <- paste("Foliose Montipora spp.")
od_mean_mospfo$sig <- c(NA,NA,NA)
od_mean_tot <- rbind(od_mean_pgwc,od_mean_pospem,od_mean_pospmd,
                        od_mean_mospem,od_mean_mospfo,od_mean_pmvc)


# od_mean_tot$DEPTH_BIN <- factor(od_mean_tot$DEPTH_BIN, levels = c("Shallow","Mid","Deep")) # Add Posthoc groupings from glms
# od_mean_tot <- od_mean[order(od_mean_tot$DEPTH_BIN),];od_mean_tot 


# Plot Partial mortality data

hist(ssss$dead)

s <- ggplot(od_mean_tot ,  aes(x=OBS_YEAR, y=Ave.od, fill=OBS_YEAR, group = Genus)) +
          geom_bar(stat = "identity") +
          scale_fill_manual(values = alpha(c("#440154FF","#22A884FF","#FDE725FF"))) +
          geom_errorbar(data = od_mean_tot,aes(ymin=Ave.od-se, ymax=Ave.od+se), width = .1) +
          geom_text(aes(x=OBS_YEAR,y=Ave.od+se+3,label=sig, group = Genus),
                    position = position_dodge()) +
          facet_wrap(~Genus,nrow = 2, labeller = label_wrap_gen(width = 18)) +
          guides(fill="none") +
           theme_bw() +
          theme(panel.grid.minor = element_blank(),
                panel.grid.major = element_blank(),
                strip.placement = "outside",
                panel.spacing = unit(0, "lines"),
                strip.text.x = element_text(size = 13, face = "bold"),
                axis.text.x = element_text(size = 13, face = "bold"),
                axis.text.y = element_text(size = 13, face = "bold"),
                axis.title.x = element_text(size=15, face = "bold"),
                axis.title.y = element_text(size=15, face = "bold")) +
          labs(x="Year",y="Average Old Dead") +
           ylim(0,39)