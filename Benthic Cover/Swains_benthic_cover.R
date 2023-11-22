rm(list = ls())

library(gdata)             # needed for drop_levels()
library(reshape)           # reshape library inclues the cast() function used below
library(splitstackshape)
library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(SimSurvey)
library(raster)
library(data.table)
library(patchwork)
library(ggdark)
library(colorRamps)
library(readr)
library(ggrepel)
library(ggnewscale)
library(ggspatial)
library(ggmap)
library(sf)
library(ggsn)
library(vegan)

source("C:/Users/jonathan.charendoff/Documents/GitHub/Benthic-Scripts/Functions/Benthic_Functions_newApp.R")
source("C:/Users/jonathan.charendoff/Documents/GitHub/fish-paste/lib/core_functions.R")
source("C:/Users/jonathan.charendoff/Documents/GitHub/fish-paste/lib/fish_team_functions.R")
source("C:/Users/jonathan.charendoff/Documents/GitHub/fish-paste/lib/Islandwide Mean&Variance Functions.R")

lu<-read.csv("T:/Benthic/Data/Lookup Tables/All_Photoquad_codes.csv")
lu$CATEGORY_CODE<-lu$TIER_1
lu$SUBCATEGORY_CODE<-lu$TIER_2b
lu$GENERA_CODE<-lu$TIER_3
lu$CODE<-lu$TIER_3

lu$CATEGORY_NAME<-lu$T1_DESC
lu$SUBCATEGORY_NAME<-lu$T2b_DESC
lu$GENERA_NAME<-lu$T3_DESC
colnames(lu)[colnames(lu)=="Cnet_SHORT_CODE"]<-"SHORT_CODE"

#FUNCTIONAL GROUP LU TABLE
new.grouping <- read.csv("C:/Users/Jonathan.Charendoff/Documents/GitHub/swains/Benthic Cover/Benthic_Cover_Lookup_Table.csv")

t1 <- read.csv("T:/Benthic/Data/REA Coral Demography & Cover/Summary Data/Site/BenthicCover_2010-2023_Tier1_SITE.csv")
t2 <- read.csv("T:/Benthic/Data/REA Coral Demography & Cover/Summary Data/Site/BenthicCover_2010-2023_Tier2b_SITE.csv")
t3 <- read.csv("T:/Benthic/Data/REA Coral Demography & Cover/Summary Data/Site/BenthicCover_2010-2023_Tier3_SITE.csv")
strata_T2 <- read.csv("T:/Benthic/Data/REA Coral Demography & Cover/Summary Data/Stratum/BenthicCover_2010-2022_Tier2_STRATA_updated.csv")

swains_t1_w <- t1[t1$ISLAND == "Swains" & t1$Oceanography == 0,-(33:36)]
swains_t2_w <- t2[t2$ISLAND == "Swains" & t2$Oceanography == 0,-(106:109)]
swains_t3_w <- t3[t3$ISLAND == "Swains" ,-(127:130)]

#wide to long format
swains_t1 <- pivot_longer(swains_t1_w, cols = c(24:ncol(swains_t1_w)), names_to = "TAXA", values_to = "PERCENT")
swains_t1$TIER <- 1

swains_t2 <- pivot_longer(swains_t2_w, cols = c(24:ncol(swains_t2_w)), names_to = "TAXA", values_to = "PERCENT")
swains_t2$TIER <- 2
test <- swains_t2 %>% 
  group_by(TAXA) %>% 
  summarise(zeros = sum(PERCENT))
zeros <- test$TAXA[test$zeros == 0]
swains_t2 <- swains_t2[-(which(swains_t2$TAXA %in% zeros)),]
swains_t2$TAXA <- recode(swains_t2$TAXA, "MONS" = "ASTS")

#NEW FUNCTIONAL GROUP LABELING
swains_t2 <- left_join(swains_t2, new.grouping)
swains_t2$New[is.na(swains_t2$New)] <- "Hard Coral"

swains_t2_pooled <- swains_t2 %>% 
  group_by(SITE, OBS_YEAR, DATE_, ANALYSIS_YEAR, new_MIN_DEPTH_M, new_MAX_DEPTH_M, LATITUDE, LONGITUDE, DEPTH_BIN, New) %>% 
  summarise(Percent = sum(PERCENT))
swains_t2_pooled_w <- pivot_wider(swains_t2_pooled, names_from = New, values_from = Percent)

swains_t3 <- pivot_longer(swains_t3_w, cols = c(24:ncol(swains_t3_w)), names_to = "TAXA", values_to = "PERCENT")
swains_t3$TIER <- 3

all_tiers <- bind_rows(swains_t1, swains_t2, swains_t3)
all_tiers <- all_tiers[-which(all_tiers$PERCENT == 0),]


#####surveyglm

library(dplyr)
library(tidyr)
library(survey)
library(multcomp)
library(emmeans)

#LOAD site-level data
swa<-filter(swains_t2_pooled,new_MAX_DEPTH_M >=3) #subset sites less than 3m?

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

site.sw$Strat_conc<-paste(site.sw$OBS_YEAR, site.sw$DEPTH_BIN,sep = "_")
site.sw$OBS_YEAR<-as.character(site.sw$OBS_YEAR)
site.sw$OBS_YEAR<-as.factor(site.sw$OBS_YEAR)
site.sw$DEPTH_BIN<-as.factor(site.sw$DEPTH_BIN)


#remove strata that have less than 2 sites
site.sw<-subset(site.sw,OBS_YEAR != "2010")

site.sw.pool <- site.sw
site.sw.pool$New <- recode_factor(site.sw.pool$New, "POSP"="Total Coral", "MOSP" = "Total Coral", "POCS" = "Total Coral",  "CCA" = "CCA", "Hard Coral" = "Total Coral", .default = "Other", "MICR" = "UPMA","UPMA" = "UPMA", "TURF" = "TURF", "EMA" = "EMA")
site.sw.pool <- site.sw.pool %>% 
  group_by(SITE,OBS_YEAR,DATE_, ANALYSIS_YEAR, new_MIN_DEPTH_M, new_MAX_DEPTH_M, LATITUDE, LONGITUDE, DEPTH_BIN, New, NH, n, sw, Strat_conc) %>% 
  summarise(Percent = sum(Percent))

#FULL MODEL
des<-svydesign(id=~1, strata=~ Strat_conc, weights=~sw,data= subset(site.sw.pool, New != "Other"))

#SUBSET SHALLOW
modR.shallow<-svyglm(Percent/100 ~ New*OBS_YEAR, design=subset(des, DEPTH_BIN == "Shallow"), family = quasibinomial(link = "logit"))
shal <- bind_cols("Shallow",car::Anova(modR.shallow, type = 3, test.statistic = "F"))
post.shal <- emmeans(modR.shallow, specs = pairwise ~ OBS_YEAR*New, adjust = "none")

#POST HOC TESTS FOR EACH FUNCTIONAL GROUP
coral.shal <- as.data.frame(post.shal[[2]]) %>% filter(stringr::str_count(contrast, 'Total Coral') == 2)
CCA.shal <- as.data.frame(post.shal[[2]]) %>% filter(stringr::str_count(contrast, 'CCA') == 2)
UPMA.shal <- as.data.frame(post.shal[[2]]) %>% filter(stringr::str_count(contrast, 'UPMA') == 2)
TURF.shal <- as.data.frame(post.shal[[2]]) %>% filter(stringr::str_count(contrast, 'TURF') == 2)
EMA.shal <- as.data.frame(post.shal[[2]]) %>% filter(stringr::str_count(contrast, 'EMA') == 2)
post.shal <- bind_rows(coral.shal, CCA.shal, UPMA.shal, TURF.shal, EMA.shal)
post.shal$p.adj <-  p.adjust(post.shal$p.value, method = "BH")
post.shal$strata <- "Shallow"

#SUBSET SHALLOW
modR.mid<-svyglm(Percent/100 ~ New*OBS_YEAR, design=subset(des, DEPTH_BIN == "Mid"), family = quasibinomial(link = "logit"))
mid <- bind_cols("Mid",car::Anova(modR.mid, type = 3, test.statistic = "F"))
post.mid <- emmeans(modR.mid, specs = pairwise ~ OBS_YEAR*New, adjust = "none")

#POST HOC TESTS FOR EACH FUNCTIONAL GROUP
coral.mid <- as.data.frame(post.mid[[2]]) %>% filter(stringr::str_count(contrast, 'Total Coral') == 2)
CCA.mid <- as.data.frame(post.mid[[2]]) %>% filter(stringr::str_count(contrast, 'CCA') == 2)
UPMA.mid <- as.data.frame(post.mid[[2]]) %>% filter(stringr::str_count(contrast, 'UPMA') == 2)
TURF.mid <- as.data.frame(post.mid[[2]]) %>% filter(stringr::str_count(contrast, 'TURF') == 2)
EMA.mid <- as.data.frame(post.mid[[2]]) %>% filter(stringr::str_count(contrast, 'EMA') == 2)
post.mid <- bind_rows(coral.mid, CCA.mid, UPMA.mid, TURF.mid, EMA.mid)
post.mid$p.adj <-  p.adjust(post.mid$p.value, method = "BH")
post.mid$strata <- "Mid"

#SUBSET SHALLOW
modR.deep<-svyglm(Percent/100 ~ New*OBS_YEAR, design=subset(des, DEPTH_BIN == "Deep"), family = quasibinomial(link = "logit"))
deep <- bind_cols("Deep", car::Anova(modR.deep, type = 3, test.statistic = "F"))
post.deep <- emmeans(modR.deep, specs = pairwise ~ OBS_YEAR*New, adjust = "none")

#POST HOC TESTS FOR EACH FUNCTIONAL GROUP
coral.deep <- as.data.frame(post.deep[[2]]) %>% filter(stringr::str_count(contrast, 'Total Coral') == 2)
CCA.deep <- as.data.frame(post.deep[[2]]) %>% filter(stringr::str_count(contrast, 'CCA') == 2)
UPMA.deep <- as.data.frame(post.deep[[2]]) %>% filter(stringr::str_count(contrast, 'UPMA') == 2)
TURF.deep <- as.data.frame(post.deep[[2]]) %>% filter(stringr::str_count(contrast, 'TURF') == 2)
EMA.deep <- as.data.frame(post.deep[[2]]) %>% filter(stringr::str_count(contrast, 'EMA') == 2)
post.deep <- bind_rows(coral.deep, CCA.deep, UPMA.deep, TURF.deep, EMA.deep)
post.deep$p.adj <-  p.adjust(post.deep$p.value, method = "BH")
post.deep$strata <- "Deep"

post.hocs <- dplyr::bind_rows(post.deep, post.mid, post.shal)
write.csv(post.hocs, "T:/Benthic/Projects/Swains 2023 Benthic Analysis/Tables/BenthicCover_posthoc.csv")
write.csv(bind_rows(shal, mid, deep), "T:/Benthic/Projects/Swains 2023 Benthic Analysis/Tables/BenthicCover_anova.csv")


##########ALL TAXA BENTHIC COVER PLOT#####
strata_T2 <- read.csv("T:/Benthic/Data/REA Coral Demography & Cover/Summary Data/Stratum/BenthicCover_2010-2022_Tier2_STRATA_updated.csv")
strata_T2 <- strata_T2[strata_T2$ISLAND == "Swains",]

#REFORMAT DATAFRAME
strata_T2_l <- pivot_longer(strata_T2, cols = c(9:176), names_to = "TAXA", values_to = "Mean")
strata_mean <- strata_T2_l %>% filter(grepl('Mean.', TAXA))
strata_mean$TAXA <- stringr::str_remove(strata_mean$TAXA, "Mean.")
strata_SE <- strata_T2_l %>% filter(grepl('SE.', TAXA))
strata_SE$TAXA <- stringr::str_remove(strata_SE$TAXA, "SE.")
colnames(strata_SE)[ncol(strata_SE)] <- "SE"
strata_T2_l <- cbind(strata_mean, SE = strata_SE$SE)
strata_T2_l$SE[is.na(strata_T2_l$SE)] <- 0

#REMOVE TAXA THAT ARE NEVER PRESENT
test <- strata_T2_l %>% 
  group_by(TAXA) %>% 
  summarise(zeros = sum(Mean))
zeros <- test$TAXA[test$zeros == 0]
strata_T2_l <- strata_T2_l[-(which(strata_T2_l$TAXA %in% zeros)),]

#RELABEL TO FUNCATIONAL GROUPS
strata_T2_l <- left_join(strata_T2_l, new.grouping)

#NEW VARIANCCE/SE WITH NEW GROUPING
strata_T2_l$variance <-(strata_T2_l$SE*sqrt(strata_T2_l$N))^2
test.pool <- strata_T2_l %>% 
  group_by(REGION,ISLAND,STRATA, ANALYSIS_YEAR, AREA_HA, N, ISLANDCODE, New) %>% 
  summarise(Mean = sum(Mean), variance = sum(variance))
strata_T2_pooled <- test.pool %>% 
  mutate(SE = sqrt(variance)/sqrt(N))

#RELABEL ANALYSIS YEAR
strata_T2_pooled$ANALYSIS_YEAR[strata_T2_pooled$ANALYSIS_YEAR == "2015-16"] <- 2015
strata_T2_pooled$ANALYSIS_YEAR <- as.Date(ISOdate(strata_T2_pooled$ANALYSIS_YEAR,1,1))
colnames(strata_T2_pooled)[8] <- "TAXA"

#SET UP TIER 1 DATA LIKE ABOVE
strata_T1 <- read.csv("T:/Benthic/Data/REA Coral Demography & Cover/Summary Data/Stratum/BenthicCover_2010-2022_Tier1_STRATA_updated.csv")
strata_T1 <- strata_T1[strata_T1$ISLAND == "Swains",]

strata_T1_l <- pivot_longer(strata_T1, cols = c(9:24), names_to = "TAXA", values_to = "Mean")
strata_mean <- strata_T1_l %>% filter(grepl('Mean.', TAXA))
strata_mean$TAXA <- stringr::str_remove(strata_mean$TAXA, "Mean.")
strata_SE <- strata_T1_l %>% filter(!grepl('Mean.', TAXA))
strata_SE$TAXA <- stringr::str_remove(strata_SE$TAXA, "SE.")
colnames(strata_SE)[ncol(strata_SE)] <- "SE"
strata_T1_l <- cbind(strata_mean, SE = strata_SE$SE)
strata_T1_l$SE[is.na(strata_T1_l$SE)] <- 0

#CALCULATE TOTAL CORAL COVER
all.coral <- strata_T1_l[strata_T1_l$TAXA == "CORAL",which(colnames(strata_T1_l) %in% colnames(strata_T2_pooled))]
all.coral$variance <- (all.coral$SE*sqrt(all.coral$N))^2
all.coral <- all.coral[,colnames(strata_T2_pooled)]

all.coral$ANALYSIS_YEAR[all.coral$ANALYSIS_YEAR == "2015-16"] <- 2015
all.coral$ANALYSIS_YEAR <- as.Date(ISOdate(all.coral$ANALYSIS_YEAR,1,1))

#ADD IN TOTAL CORAL TO TIER 2 DATA (TO HAVE BOTH FOCAL TAXA AND TOTAL CORAL IN THE SAME DF)
strata_T2_pooled <- bind_rows(strata_T2_pooled, all.coral)
strata_T2_pooled$TAXA[strata_T2_pooled$TAXA == "CORAL"] <- "Total Coral"
taxalist <- c("Total Coral", "POSP", "MOSP", "POCS", 
              "CCA","UPMA" , "TURF","EMA", 
              "Hard Coral", "MICR", "TUN", "Other_Invert", "Non-Biological")

#RELABEL FACTOR FOR PLOTS
strata_T2_pooled$TAXA <- factor(strata_T2_pooled$TAXA, levels = taxalist)
strata_T2_pooled$STRATA <- factor(strata_T2_pooled$STRATA, levels = c("FS", "FM", "FD"))
strata_T2_pooled$STRATA <- recode_factor(strata_T2_pooled$STRATA, "FS"="Shallow", "FM" = "Mid", "FD" = "Deep")
strata_T2_pooled$TAXA <- recode_factor(strata_T2_pooled$TAXA, "POSP"="Porites sp.", "MOSP" = "Montipora sp.", "POCS" = "Pocillopora sp.", "TUN" = "Tunicate", "Other_Invert" = "Other Invert.")
taxacolors <- c("black",  
                "hotpink1", "palegreen4", "limegreen","orange2", 
                "coral", "lightgreen", "lightskyblue1", "lightcyan3", "azure4")

strata_T2_pooled$Year <- as.character(lubridate::year(strata_T2_pooled$ANALYSIS_YEAR))


lines <- ggplot(subset(strata_T2_pooled, TAXA %in% c("Total Coral","CCA","UPMA", "TURF", "EMA")), aes(x = ANALYSIS_YEAR, y = Mean, color = TAXA)) + 
  geom_point(size = 4)+ 
  geom_line()+
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 1)+
  facet_wrap(~STRATA, scales = "free_y", ncol = 1)+
  scale_x_date(date_labels = "%Y", 
               date_breaks = "1 year", 
               limits = as.Date(c("2009-06-01", "2023-12-01")),
               expand = c(0,0)
  )+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        strip.background = element_blank())+
  labs(x = "Year",
       y = "Percent Cover",
       color = "Functional Group")+
  scale_color_manual(values = taxacolors)+
  geom_vline(xintercept = as.Date("2016-01-01"), color = "red", linetype = 2)


ggsave("T:/Benthic/Projects/Swains 2023 Benthic Analysis/plots/benthic_change.png", lines, height = 8, width = 11)

#####STR DATA#####
temps_2015 <- read.csv("C:/Users/Jonathan.Charendoff/Desktop/ESD_NCRMP_Temperature_2015_AMSM.csv")
temps_2018 <- read.csv("C:/Users/Jonathan.Charendoff/Desktop/ESD_NCRMP_Temperature_2018_AMSM.csv")
swains.temps.2015 <- temps_2015[temps_2015$LOCATION == "Swains",]
swains.temps.2018 <- temps_2018[temps_2018$LOCATION == "Swains",]
swains.temps.2018$INSTRUMENTSN <- as.character(swains.temps.2018$INSTRUMENTSN)
swains.temps <- bind_rows(swains.temps.2015, swains.temps.2018)
swains.temps$DEPTH_BIN[swains.temps$OCC_SITEID == "OCC-SWA-002"] <- "Deep"
swains.temps$DEPTH_BIN[swains.temps$OCC_SITEID == "OCC-SWA-003"] <- "Mid"
swains.temps$DEPTH_BIN[swains.temps$OCC_SITEID == "OCC-SWA-004"] <- "Shallow"
swains.temps$DEPTH_BIN[swains.temps$OCC_SITEID == "OCC-SWA-006"] <- "Deep"
swains.temps$DEPTH_BIN[swains.temps$OCC_SITEID == "OCC-SWA-007"] <- "Mid"
swains.temps$DEPTH_BIN[swains.temps$OCC_SITEID == "OCC-SWA-008"] <- "Shallow"
swains.temps$DEPTH_BIN[swains.temps$OCC_SITEID == "OCC-SWA-009"] <- "Shallow"
swains.temps$DEPTH_BIN <- factor(swains.temps$DEPTH_BIN, levels = c("Shallow", "Mid", "Deep"))


swains.temps.monthly <- swains.temps %>% 
  group_by(DEPTH_BIN, YEAR, MONTH) %>% 
  summarise(TEMP_C = mean(TEMP_C))
swains.temps.monthly$DATE <- as.Date(paste(swains.temps.monthly$YEAR, swains.temps.monthly$MONTH, "1",sep = "-"), format = "%Y-%m-%d")

colors <- c("#AA4499" ,"#44BB99", "#FF8C00")

STR_temps <- ggplot(swains.temps.monthly[swains.temps.monthly$YEAR %in% c(2015,2016),], aes(x = DATE, y = TEMP_C, color = DEPTH_BIN))+
  geom_line()+
  scale_color_manual(values = colors)+
  theme_bw()+
  labs(x = "Date",
       y = expression("Temperature " ( degree~C)),
       color = "Depth Bin")+
  theme(legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16))

ggsave("T:/Benthic/Projects/Swains 2023 Benthic Analysis/plots/STR_TEMPS.png", STR_temps, height = 8, width = 11)


######Calcifiers#####
site.sw.calc <- site.sw
site.sw.calc$New <- recode_factor(site.sw.calc$New, "POSP"="CALC", "MOSP" = "CALC", "POCS" = "CALC",  "CCA" = "CALC", "Hard Coral" = "CALC", .default = "NC")
site.sw.calc <- site.sw.calc %>% 
  group_by(SITE,OBS_YEAR,DATE_, ANALYSIS_YEAR, new_MIN_DEPTH_M, new_MAX_DEPTH_M, LATITUDE, LONGITUDE, DEPTH_BIN, New, NH, n, sw, Strat_conc) %>% 
  summarise(Percent = sum(Percent)/100)

site.sw.calc.w <- pivot_wider(site.sw.calc, names_from = New, values_from = Percent)

site.sw.calc <- droplevels(site.sw.calc %>% filter(New == "CALC")) # 

des.calc<-svydesign(id=~1, strata=~ Strat_conc, weights=~sw,data=site.sw.calc)

hist(car::logit(site.sw.calc$Percent))
with(site.sw.calc, tapply((car::logit(Percent)), OBS_YEAR, shapiro.test)) #confirmed normal distribution among factor levels within OBS_YEAR
bartlett.test((car::logit(Percent)) ~ OBS_YEAR, site.sw.calc) #confirmed homogeneity of variance among factor levels within OBS_YEAR
with(site.sw.calc, tapply((car::logit(Percent)), DEPTH_BIN, shapiro.test)) #confirmed normal distribution for DEPTH_BIN
bartlett.test((car::logit(Percent)) ~ DEPTH_BIN, site.sw.calc) #confirmed homogeneity of variance for DEPTH_BIN

modR.calc<-svyglm(Percent ~ OBS_YEAR*DEPTH_BIN, design=des.calc, family = quasibinomial(link = "logit"))
summary(modR.calc)
calc.anova <- car::Anova(modR.calc, type = 3, test.statistic = "F")
post.calc <- summary(glht(modR.calc, mcp(DEPTH_BIN="Tukey")))
write.csv(calc.anova, "T:/Benthic/Projects/Swains 2023 Benthic Analysis/Tables/CalcCover.csv")

##calcifiers vs non calcifiers
calcified <- strata_T1_l
calcified$CALC <- recode_factor(calcified$TAXA, "CORAL"="Calcifier", "CCA" = "Calcifier", .default = "Non-Calcifier")
calcified$variance <-(calcified$SE*sqrt(calcified$N))^2
calc.pool <- calcified %>% 
  group_by(REGION,ISLAND,STRATA, ANALYSIS_YEAR, AREA_HA, N, ISLANDCODE, CALC) %>% 
  summarise(Mean = sum(Mean), variance = sum(variance))
calc_pooled <- calc.pool %>% 
  mutate(SE = sqrt(variance)/sqrt(N))
calc_pooled$ANALYSIS_YEAR[calc_pooled$ANALYSIS_YEAR == "2015-16"] <- 2015
calc_pooled$ANALYSIS_YEAR <- as.Date(ISOdate(calc_pooled$ANALYSIS_YEAR,1,1))
calc_pooled$STRATA <- factor(calc_pooled$STRATA, levels = c("FS", "FM", "FD"))
calc_pooled$STRATA <- recode_factor(calc_pooled$STRATA, "FS"="Shallow", "FM" = "Mid", "FD" = "Deep")

calc <- ggplot(calc_pooled, aes(x = ANALYSIS_YEAR, y = Mean, color = CALC)) + 
  geom_point(size = 4)+ 
  geom_line()+
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 1)+
  facet_wrap(~STRATA, scales = "free_y", ncol = 1)+
  scale_x_date(date_labels = "%Y", 
               date_breaks = "1 year", 
               limits = as.Date(c("2009-06-01", "2023-12-01")),
               expand = c(0,0)
  )+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        strip.background = element_blank())+
  labs(x = "Year",
       y = "Percent Cover",
       color = "Functional Group")+
  scale_color_manual(values = c("hotpink1", "palegreen4"))+
  geom_vline(xintercept = as.Date("2016-01-01"), color = "red", linetype = 2)

ggsave("T:/Benthic/Projects/Swains 2023 Benthic Analysis/plots/calcifier_change.png", calc, height = 8, width = 11)




######Focal Taxa######
#Subset data
site.sw.focal <- site.sw %>% filter(New == "POCS" |New == "POSP"| New == "MOSP") %>% droplevels()
site.sw.focal$New <- as.factor(site.sw.focal$New)

des.focal<-svydesign(id=~1, strata=~ Strat_conc, weights=~sw,data=site.sw.focal)


hist(car::logit(site.sw.posp$Percent/100))
with(site.sw.pocs, tapply((log(Percent+1)), OBS_YEAR, shapiro.test)) #confirmed normal distribution among factor levels within OBS_YEAR
bartlett.test((log(Percent+1)) ~ OBS_YEAR, site.sw.pocs) #confirmed homogeneity of variance among factor levels within OBS_YEAR
with(site.sw.pocs, tapply((log(Percent+1)), DEPTH_BIN, shapiro.test)) #confirmed normal distribution for DEPTH_BIN
bartlett.test((log(Percent+1)) ~ DEPTH_BIN, site.sw.pocs) #confirmed homogeneity of variance for DEPTH_BIN
#log transform pocs


with(site.sw.mosp, tapply((sqrt(Percent)), OBS_YEAR, shapiro.test)) #confirmed normal distribution among factor levels within OBS_YEAR
bartlett.test((sqrt(Percent)) ~ OBS_YEAR, site.sw.mosp) #confirmed homogeneity of variance among factor levels within OBS_YEAR
with(site.sw.mosp, tapply((sqrt(Percent)), DEPTH_BIN, shapiro.test)) #confirmed normal distribution for DEPTH_BIN
bartlett.test((sqrt(Percent)) ~ DEPTH_BIN, site.sw.mosp) #confirmed homogeneity of variance for DEPTH_BIN
#sqrt transform mosp

with(site.sw.posp, tapply((log(Percent+1)), OBS_YEAR, shapiro.test)) #confirmed normal distribution among factor levels within OBS_YEAR
bartlett.test((log(Percent+1)) ~ OBS_YEAR, site.sw.posp) #confirmed homogeneity of variance among factor levels within OBS_YEAR
with(site.sw.posp, tapply((log(Percent+1)), DEPTH_BIN, shapiro.test)) #confirmed normal distribution for DEPTH_BIN
bartlett.test((log(Percent+1)) ~ DEPTH_BIN, site.sw.posp) #confirmed homogeneity of variance for DEPTH_BIN
#log transform posp

modR.coral.shal<-svyglm(Percent/100 ~ New*OBS_YEAR, design=subset(des.focal, DEPTH_BIN == "Shallow"), family = quasibinomial(link = "logit"))
modR.coral.mid<-svyglm(Percent/100 ~ New*OBS_YEAR, design=subset(des.focal, DEPTH_BIN == "Mid"), family = quasibinomial(link = "logit")) # Depth and year:depth significant
modR.coral.deep<-svyglm(Percent/100 ~ New*OBS_YEAR, design=subset(des.focal, DEPTH_BIN == "Deep"), family = quasibinomial(link = "logit")) # Depth significant

coral.shal <- bind_cols("shallow",car::Anova(modR.coral.shal, type = 3, test.statistic = "F")) #only depth
coral.mid <- bind_cols("mid",car::Anova(modR.coral.mid, type = 3, test.statistic = "F"))#only depth
coral.deep <- bind_cols("deep",car::Anova(modR.coral.deep, type = 3, test.statistic = "F")) #only depth
write.csv(bind_rows(coral.shal, coral.mid, coral.deep), "T:/Benthic/Projects/Swains 2023 Benthic Analysis/Tables/BenthicCoralCover_anova.csv")

post.coral.shal <- emmeans(modR.coral.shal, specs = pairwise ~ OBS_YEAR*New, adjust = "none")
post.coral.mid <- emmeans(modR.coral.mid, specs = pairwise ~ OBS_YEAR*New, adjust = "none")
post.coral.deep <- emmeans(modR.coral.deep, specs = pairwise ~ OBS_YEAR*New, adjust = "none")

POCS.shal <- as.data.frame(post.coral.shal[[2]]) %>% filter(stringr::str_count(contrast, 'POCS') == 2)
MOSP.shal <- as.data.frame(post.coral.shal[[2]]) %>% filter(stringr::str_count(contrast, 'MOSP') == 2)
POSP.shal <- as.data.frame(post.coral.shal[[2]]) %>% filter(stringr::str_count(contrast, 'POSP') == 2)

post.coral.shal <- bind_rows(POCS.shal, MOSP.shal, POSP.shal)
post.coral.shal$p.adj <-  p.adjust(post.coral.shal$p.value)
post.coral.shal$strata <- "Shallow"

POCS.mid <- as.data.frame(post.coral.mid[[2]]) %>% filter(stringr::str_count(contrast, 'POCS') == 2)
MOSP.mid <- as.data.frame(post.coral.mid[[2]]) %>% filter(stringr::str_count(contrast, 'MOSP') == 2)
POSP.mid <- as.data.frame(post.coral.mid[[2]]) %>% filter(stringr::str_count(contrast, 'POSP') == 2)

post.coral.mid <- bind_rows(POCS.mid, MOSP.mid, POSP.mid)
post.coral.mid$p.adj <-  p.adjust(post.coral.mid$p.value)
post.coral.mid$strata <- "Mid"

post.hocs.coral <- bind_rows(post.coral.mid, post.coral.shal)
write.csv(post.hocs.coral, "T:/Benthic/Projects/Swains 2023 Benthic Analysis/Tables/BenthicCoralCover_posthoc.csv")

corals <- ggplot(subset(strata_T2_pooled, TAXA %in% c("Porites sp.","Montipora sp.","Pocillopora sp.")), aes(x = ANALYSIS_YEAR, y = Mean, color = TAXA)) + 
  geom_point(size = 4)+ 
  geom_line()+
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 1)+
  facet_wrap(~STRATA, scales = "free_y", ncol = 1)+
  scale_x_date(date_labels = "%Y", 
               date_breaks = "1 year", 
               limits = as.Date(c("2009-06-01", "2023-12-01")),
               expand = c(0,0)
  )+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        strip.background = element_blank())+
  labs(x = "Year",
       y = "Percent Cover",
       color = "Depth Strata")+
  geom_vline(xintercept = as.Date("2016-01-01"), color = "red", linetype = 2)+
  scale_color_manual(values = c("#00cdd0","#ff8a01","#006400"))

ggsave("T:/Benthic/Projects/Swains 2023 Benthic Analysis/plots/coral_change.png", corals, height = 8, width = 11)

#########map with sites######

shape <- st_read(dsn = "T:/Common/Maps/Island", layer = "islands")
shape <- subset(shape, ISLAND == "SWAINS")

all.sites$OBS_YEAR <- as.factor(all.sites$OBS_YEAR)
swains_t2_w$Calcifiers <- swains_t1_w$CCA + swains_t1_w$CORAL +swains_t1_w$HAL

maps<- ggplot() + 
  geom_sf(data = shape, inherit.aes = F) + 
  geom_spatial_point(data = subset(swains_t2_pooled, New %in% c("MOSP", "POCS", "POSP")), aes(x = LONGITUDE, y = LATITUDE, color = Percent))+
  facet_grid(New~OBS_YEAR)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.text.x = element_text(angle = 90))+
  scale_color_viridis_c(direction = -1)+
  labs(
    x = NULL, 
    y = NULL, 
    title = "All Sites")

ggsave("T:/Benthic/Projects/Swains 2023 Benthic Analysis/plots/Calc_Cover_by_DepthxYear.png", maps, height = 8, width = 11)
