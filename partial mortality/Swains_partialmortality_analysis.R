#Data do not include FUSP

# Using R version 4.1.0 (2021-05-18)

rm(list=ls())

library(dplyr)
library(tidyr)
library(survey)
library(multcomp)
library(ggplot2)
library(emmeans)

#LOAD site-level data
site.data.gen2<-read.csv("C:/Users/Corinne.Amir/Documents/GitHub/swains/Data/Swains_sitedata_GENUS.csv") # For genus results
swam <- read.csv("C:/Users/Corinne.Amir/Documents/GitHub/swains/Data/Swains_sitedata_TAXONCODE_MORPH.csv") # For simper results


#### Genus level ####
swa<-dplyr::filter(site.data.gen2,new_MAX_DEPTH_M >=3) #subset sites less than 3m?

nrow(swa)

#read in sector-area file
sectors<-read.csv("C:/Users/Corinne.Amir/Documents/GitHub/swains/Data/Sectors-Strata-Areas.csv", stringsAsFactors=FALSE)
swa_sa<-dplyr::filter(sectors,ISLAND=="Swains")
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

nrow(site.sw)



#Create contactenated Strata variable 
# site.sw$Strat_conc<-paste(site.sw$OBS_YEAR, site.sw$DB_RZ,sep = "_")
site.sw$Strat_conc<-paste(site.sw$OBS_YEAR, site.sw$DEPTH_BIN,sep = "_")
site.sw$OBS_YEAR<-as.factor(site.sw$OBS_YEAR)


#Subset data (selection based on most common genera)
site.sw.pocs <- site.sw %>% filter(GENUS_CODE == "POCS") # above code isn't 50 so here I'm looking at just total coral
site.sw.mosp <- site.sw %>% filter(GENUS_CODE == "MOSP") # 
site.sw.posp <- site.sw %>% filter(GENUS_CODE == "POSP") # 

#remove NA values from POSP df (can't do this in the svyby function)
site.sw.posp<-site.sw.posp[complete.cases(site.sw.posp), ]

site.sw <- site.sw %>% filter(GENUS_CODE == "SSSS")

#remove strata that have less than 2 sites
site.sw<-subset(site.sw,n>1)
summary(site.sw$n)
length(unique(site.sw$SITE))

site.sw.pocs<-subset(site.sw.pocs,n>1)
site.sw.mosp<-subset(site.sw.mosp,n>1)
site.sw.posp<-subset(site.sw.posp,n>1)
site.sw.posp<-subset(site.sw.posp, Strat_conc !="2018_Mid")

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
depth_mean<-svyby(~Ave.od,~OBS_YEAR+DEPTH_BIN,des,svymean)
depth_mean_pocs<-svyby(~Ave.od,~OBS_YEAR+DEPTH_BIN,des.pocs,svymean)
depth_mean_mosp<-svyby(~Ave.od,~OBS_YEAR+DEPTH_BIN,des.mosp,svymean)
depth_mean_posp<-svyby(~Ave.od,~OBS_YEAR+DEPTH_BIN,des.posp,svymean)
# depth_mean<-svyby(~Ave.od,~OBS_YEAR+DEPTH_BIN,des,svymean)

#Testing for normality and equal variance
with(site.sw, tapply((Ave.od), OBS_YEAR, shapiro.test)) #confirmed normal distribution among factor levels within OBS_YEAR
bartlett.test(Ave.od ~ OBS_YEAR, site.sw) #confirmed homogeneity of variance among factor levels within OBS_YEAR
with(site.sw, tapply((Ave.od), DEPTH_BIN, shapiro.test)) #confirmed normal distribution for DEPTH_BIN
bartlett.test(Ave.od ~ DEPTH_BIN, site.sw) #confirmed homogeneity of variance for DEPTH_BIN
#SSSS OD - no problems

with(site.sw.pocs, tapply((sqrt(Ave.od)), OBS_YEAR, shapiro.test)) #confirmed normal distribution among factor levels within OBS_YEAR
bartlett.test((sqrt(Ave.od)) ~ OBS_YEAR, site.sw.pocs) #confirmed homogeneity of variance among factor levels within OBS_YEAR
with(site.sw.pocs, tapply((sqrt(Ave.od)), DEPTH_BIN, shapiro.test)) #confirmed normal distribution for DEPTH_BIN
bartlett.test((sqrt(Ave.od)) ~ DEPTH_BIN, site.sw.pocs) #confirmed homogeneity of variance for DEPTH_BIN
#need so sqrt transform pocs OD


with(site.sw.mosp, tapply((Ave.od), OBS_YEAR, shapiro.test)) #confirmed normal distribution among factor levels within OBS_YEAR
bartlett.test((Ave.od) ~ OBS_YEAR, site.sw.mosp) #confirmed homogeneity of variance among factor levels within OBS_YEAR
with(site.sw.mosp, tapply((Ave.od), DEPTH_BIN, shapiro.test)) #confirmed normal distribution for DEPTH_BIN
bartlett.test((Ave.od) ~ DEPTH_BIN, site.sw.mosp) #confirmed homogeneity of variance for DEPTH_BIN
#no transformation needed

with(site.sw.posp, tapply((Ave.od), OBS_YEAR, shapiro.test)) #confirmed normal distribution among factor levels within OBS_YEAR
bartlett.test((Ave.od) ~ OBS_YEAR, site.sw.posp) #confirmed homogeneity of variance among factor levels within OBS_YEAR
with(site.sw.posp, tapply((Ave.od), DEPTH_BIN, shapiro.test)) #confirmed normal distribution for DEPTH_BIN
bartlett.test((Ave.od) ~ DEPTH_BIN, site.sw.posp) #confirmed homogeneity of variance for DEPTH_BIN
#no transformation needed


#  Test fixed effects of region and year
modR<-svyglm(Ave.od ~ OBS_YEAR*DEPTH_BIN, design=des) # all significant
modR.pocs<-svyglm((sqrt(Ave.od)) ~ OBS_YEAR*DEPTH_BIN, design=des.pocs) # Year and depth significant
modR.mosp<-svyglm(Ave.od ~ OBS_YEAR*DEPTH_BIN, design=des.mosp) # Depth and year:depth significant
modR.posp<-svyglm(Ave.od ~ OBS_YEAR*DEPTH_BIN, design=des.posp) # Depth significant

summary(modR); summary(modR.pocs); summary(modR.mosp); summary(modR.posp)
car::Anova(modR, type = 3, test.statistic = "F") #nothing significant
car::Anova(modR.pocs, type = 3, test.statistic = "F") #only year signficant
car::Anova(modR.mosp, type = 3, test.statistic = "F") #nothing sig
car::Anova(modR.posp, type = 3, test.statistic = "F") #all factors signficant


# Post-hoc
# Year
y <- svyglm((sqrt(Ave.od)) ~ OBS_YEAR*DEPTH_BIN, design=des.pocs) # All corals
emmeans(y, pairwise ~ OBS_YEAR) #od increased from 2015 and 2018

y <- svyglm(Ave.od ~ OBS_YEAR, design = subset(des.posp,DEPTH_BIN=="Shallow")) # POCS
emmeans(y, pairwise ~ OBS_YEAR) #od increased from 2015 and 2018
            
y <- svyglm(Ave.od ~ OBS_YEAR, design = subset(des.posp,DEPTH_BIN=="Mid")) # POCS
emmeans(y, pairwise ~ OBS_YEAR) #od increased from 2015 and 2018 then decreased from 2018-2023

y <- svyglm(Ave.od ~ OBS_YEAR, design = subset(des.posp,DEPTH_BIN=="Deep")) # POCS
emmeans(y, pairwise ~ OBS_YEAR) #no changes at deep sites

summary(glht(modR.pocs, mcp(OBS_YEAR="Tukey"))) # 2023 sig different from 2015
summary(glht(modR.mosp, mcp(OBS_YEAR="Tukey"))) # 2023 and 2018 sig different than 2015
summary(glht(modR.posp, mcp(OBS_YEAR="Tukey"))) # 2023 and 2018 sig different than 2015


# Put all depth_mean dataframes together for plotting
#depth_mean$Genus <- paste("Total")
depth_mean_pocs$Genus <- paste("Pocillopora")
depth_mean_mosp$Genus <- paste("Montipora")
#depth_mean_posp$Genus <- paste("Porites")
depth_mean_tot <- rbind(depth_mean_pocs, depth_mean_mosp)

# Plot Data
s <- ggplot(depth_mean_tot %>% filter(DEPTH_BIN == "Shallow"),
            aes(x = OBS_YEAR, y = Ave.od, group = Genus, fill = OBS_YEAR)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = alpha(c("#440154FF","#22A884FF","#FDE725FF"))) +
  geom_errorbar(data = depth_mean_tot %>% filter(DEPTH_BIN == "Shallow"),
                aes(ymin = Ave.od-se, ymax = Ave.od+se),
                width = .2) +
  facet_wrap(~Genus, nrow = 1) +
  guides(fill="none") +
  
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
  # geom_text(aes(x=OBS_YEAR,y=Ave.od+se,label=c("","ns","","","ns","","a","b","b")),
  #           position = position_dodge(),
  #           vjust = -0.5) 

m <- ggplot(depth_mean_tot %>% filter(DEPTH_BIN == "Mid"), 
            aes(x = OBS_YEAR, y = Ave.od, group = Genus, fill = OBS_YEAR)) +
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = alpha(c("#440154FF","#22A884FF","#FDE725FF"))) +
  geom_errorbar(data = depth_mean_tot %>% filter(DEPTH_BIN == "Mid"),
                aes(ymin = Ave.od-se, ymax = Ave.od+se),
                width = .1) +
  facet_wrap(~Genus, nrow = 1) +
  guides(fill="none") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
  # geom_text(aes(x=OBS_YEAR,y=Ave.od+se,label=c("","ns","","","ns","","a","a","b")),
  #           position = position_dodge(),
  #           vjust = -0.5) 

d <- ggplot(depth_mean_tot %>% filter(DEPTH_BIN == "Deep"), 
            aes(x = OBS_YEAR, y = Ave.od, group = Genus, fill = OBS_YEAR)) +
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = alpha(c("#440154FF","#22A884FF","#FDE725FF"))) +
  geom_errorbar(data = depth_mean_tot %>% filter(DEPTH_BIN == "Deep"),
                aes(ymin = Ave.od-se, ymax = Ave.od+se),
                width = .1) +
  facet_wrap(~Genus, nrow = 1) +
  guides(fill="none") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
  # geom_text(aes(x=OBS_YEAR,y=Ave.od+se,label=c("a","b","b","","ns","","a","b","a")),
  #           position = position_dodge(),
  #           vjust = -0.5)   

tot <- ggplot(depth_mean, 
              aes(x = OBS_YEAR, y = Ave.od, fill = OBS_YEAR)) +
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = alpha(c("#440154FF","#22A884FF","#FDE725FF"))) +
  geom_errorbar(data = depth_mean,
                aes(ymin = Ave.od-se, ymax = Ave.od+se),
                width = .1) +
  facet_wrap(~DEPTH_BIN, nrow = 1) +
  guides(fill="none") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) 
ytitle="Mean Old Dead"
xtitle= "Year"

grid.arrange(s+ggtitle("Shallow"), m + ggtitle("Mid"), d+ ggtitle("Deep"), 
             left = ytitle,
             bottom = xtitle,
             nrow = 3)



#### SIMPER taxa level ####
# Merge demography data with new NH values pooled across the 2 swains sectors
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


s <- ggplot(od_mean_tot ,  aes(x=OBS_YEAR, y=Ave.od, fill=OBS_YEAR, group = Genus)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = alpha(c("orchid4","aquamarine3","goldenrod3"))) +
  geom_errorbar(data = od_mean_tot,aes(ymin=Ave.od-se, ymax=Ave.od+se), width = .3) +
  geom_text(aes(x=OBS_YEAR,y=Ave.od+se+3,label=sig, group = Genus),
            position = position_dodge()) +
  facet_wrap(~Genus,nrow = 2, labeller = label_wrap_gen(width = 18)) +
  guides(fill="none") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        strip.placement = "outside",
        panel.spacing = unit(0, "lines"),
        strip.text.x = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size=13),
        axis.title.y = element_text(size=13)) +
  labs(x="Year",y="Average Old Dead %") +
  ylim(0,39)

ggsave(plot=s,file="C:/Users/Corinne.Amir/Documents/Analysis/FigureS5.jpg",width=10,height=8)
ggsave(plot=s,file="T:/Benthic/Projects/Swains 2023 Benthic Analysis/Plots/FigureS5.jpg",width=10,height=8)
