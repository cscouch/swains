rm(list = ls())

library(gdata)             # needed for drop_levels()
library(reshape)           # reshape library inclues the cast() function used below
library(splitstackshape)
library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(survey)
library(multcomp)
library(emmeans)

dir = Sys.info()[7]
setwd(paste0("C:/Users/", dir, "/Documents/GitHub/swains/"))

#####surveyglm
#LOAD site-level data
swains_t2_pooled <- read.csv("Data/Benthic_Cover_Ready.csv")
swa<-filter(swains_t2_pooled,new_MAX_DEPTH_M >=3) #subset sites deeper than 3m

nrow(swa)

#read in sector-area file
sectors<-read.csv("Data/Sectors-Strata-Areas.csv", stringsAsFactors=FALSE)
swa_sa<-filter(sectors,ISLAND=="Swains") #subset to swains only
swa_sa$DEPTH_BIN<-as.factor(swa_sa$DEPTH_BIN)

NH <- swa_sa %>%
  group_by(SEC_NAME, DEPTH_BIN)%>%
  dplyr::summarize(unique = NH)%>%
  group_by(DEPTH_BIN)%>%
  dplyr::summarise(NH = sum(unique)) #concatenate swains open and sanctuary sectors bc we don't care about differences between these areas

swa<-left_join(swa,NH) #merge percent cover data with new NH values pooled across the 2 swains sectors


#Calculate survey weights (inverse proportion weighting)
w.df<-swa %>%
  group_by(OBS_YEAR,DEPTH_BIN,NH) %>%
  dplyr::summarise(n = length(unique(SITE)))

w.df$sw<-w.df$NH/w.df$n #calculate survey weights for each site

site.sw<-left_join(swa,w.df) #merge weights with site-level data
head(site.sw)

#site.sw$Strat_conc<-paste(site.sw$OBS_YEAR, site.sw$DEPTH_BIN,sep = "_")
site.sw$OBS_YEAR<-as.character(site.sw$OBS_YEAR)
site.sw$OBS_YEAR<-as.factor(site.sw$OBS_YEAR)
site.sw$DEPTH_BIN<-as.factor(site.sw$DEPTH_BIN)

#create strata levels for survey design
site.sw.pool <- site.sw
site.sw.pool$Strat_conc<-paste(site.sw.pool$OBS_YEAR, site.sw.pool$DEPTH_BIN, site.sw.pool$New, sep = "_")#create strata design including taxa so we can survey lonely.psu correctly

site.sw.pool.shal <- droplevels(filter(site.sw.pool, New != "Other", DEPTH_BIN == "Shallow"))
site.sw.pool.mid <- droplevels(filter(site.sw.pool, New != "Other", DEPTH_BIN == "Mid"))
site.sw.pool.deep <- droplevels(filter(site.sw.pool, New != "Other", DEPTH_BIN == "Deep"))

options(survey.lonely.psu = "adjust") #keep in strata with only 1 site

#create survey design

###CREATE DESIGNS
des.shal<-svydesign(id=~1, strata=~ Strat_conc, weights=~sw,data= site.sw.pool.shal) #create survey design including all depth strata, remove the random taxa group bc we aren't looking for differences between groups
des.mid<-svydesign(id=~1, strata=~ Strat_conc, weights=~sw,data= site.sw.pool.mid)
des.deep<-svydesign(id=~1, strata=~ Strat_conc, weights=~sw,data= site.sw.pool.deep)

#SUBSET SHALLOW
modR.shallow<-svyglm(Percent/100 ~ New*OBS_YEAR, design=des.shal, family = quasibinomial(link = "logit")) #run survey model
shal <- bind_cols("Shallow",car::Anova(modR.shallow, type = 3, test.statistic = "F")) #create anova table and add in the strata label
post.shal <- emmeans(modR.shallow, specs = pairwise ~ OBS_YEAR*New, adjust = "none") #Calculate pairwise comparisons, without correction

#Single Factor tests -- all are significant
#shal.coral <- car::Anova(svyglm(Percent/100 ~ OBS_YEAR, design=subset(des, DEPTH_BIN == "Shallow" & New == "Total Coral"), family = quasibinomial(link = "logit")), type = 3, test.statistic = "F")#p = 0.003
#shal.CCA <- car::Anova(svyglm(Percent/100 ~ OBS_YEAR, design=subset(des, DEPTH_BIN == "Shallow" & New == "CCA"), family = quasibinomial(link = "logit")), type = 3, test.statistic = "F")#p = 0.006
#shal.UPMA <- car::Anova(svyglm(Percent/100 ~ OBS_YEAR, design=subset(des, DEPTH_BIN == "Shallow" & New == "UPMA"), family = quasibinomial(link = "logit")), type = 3, test.statistic = "F")#p = 0.003
#shal.TURF <- car::Anova(svyglm(Percent/100 ~ OBS_YEAR, design=subset(des, DEPTH_BIN == "Shallow" & New == "TURF"), family = quasibinomial(link = "logit")), type = 3, test.statistic = "F")#p = 0.0065
#shal.EMA <- car::Anova(svyglm(Percent/100 ~ OBS_YEAR, design=subset(des, DEPTH_BIN == "Shallow" & New == "EMA"), family = quasibinomial(link = "logit")), type = 3, test.statistic = "F")#p = 0.04

#POST HOC TESTS FOR EACH FUNCTIONAL GROUP ONLY FOR COMPARISONS WE CARE ABOUT (WITHOUT FUNCTIONAL GROUP OVER TIME, NOT BETWEEN FUNCTIONAL GROUPS)
coral.shal <- as.data.frame(post.shal[[2]]) %>% filter(stringr::str_count(contrast, 'Total Coral') == 2)
CCA.shal <- as.data.frame(post.shal[[2]]) %>% filter(stringr::str_count(contrast, 'CCA') == 2)
UPMA.shal <- as.data.frame(post.shal[[2]]) %>% filter(stringr::str_count(contrast, 'UPMA') == 2)
TURF.shal <- as.data.frame(post.shal[[2]]) %>% filter(stringr::str_count(contrast, 'TURF') == 2)
EMA.shal <- as.data.frame(post.shal[[2]]) %>% filter(stringr::str_count(contrast, 'EMA') == 2)

#apply test correction to each functional group
coral.shal$p.adj <-  p.adjust(coral.shal$p.value, method = "BH") 
CCA.shal$p.adj <-  p.adjust(CCA.shal$p.value, method = "BH") 
UPMA.shal$p.adj <-  p.adjust(UPMA.shal$p.value, method = "BH") 
TURF.shal$p.adj <-  p.adjust(TURF.shal$p.value, method = "BH") 
EMA.shal$p.adj <-  p.adjust(EMA.shal$p.value, method = "BH") 

post.shal <- bind_rows(coral.shal, CCA.shal, UPMA.shal, TURF.shal, EMA.shal) #put the funtional groups back together
post.shal$strata <- "Shallow" #add in strata label for further grouping

#SUBSET MID 
modR.mid<-svyglm(Percent/100 ~ New*OBS_YEAR, design=des.mid, family = quasibinomial(link = "logit"))
mid <- bind_cols("Mid",car::Anova(modR.mid, type = 3, test.statistic = "F"))
post.mid <- emmeans(modR.mid, specs = pairwise ~ OBS_YEAR*New, adjust = "none")

#Single Factor tests
#mid.coral <- car::Anova(svyglm(Percent/100 ~ OBS_YEAR, design=subset(des, DEPTH_BIN == "Mid" & New == "Total Coral"), family = quasibinomial(link = "logit")), type = 3, test.statistic = "F")#p = 0.048
#mid.CCA <- car::Anova(svyglm(Percent/100 ~ OBS_YEAR, design=subset(des, DEPTH_BIN == "Mid" & New == "CCA"), family = quasibinomial(link = "logit")), type = 3, test.statistic = "F")#NS p = 0.12
#mid.UPMA <- car::Anova(svyglm(Percent/100 ~ OBS_YEAR, design=subset(des, DEPTH_BIN == "Mid" & New == "UPMA"), family = quasibinomial(link = "logit")), type = 3, test.statistic = "F")#p = 0.0001
#mid.TURF <- car::Anova(svyglm(Percent/100 ~ OBS_YEAR, design=subset(des, DEPTH_BIN == "Mid" & New == "TURF"), family = quasibinomial(link = "logit")), type = 3, test.statistic = "F")#p = 0.003
#mid.EMA <- car::Anova(svyglm(Percent/100 ~ OBS_YEAR, design=subset(des, DEPTH_BIN == "Mid" & New == "EMA"), family = quasibinomial(link = "logit")), type = 3, test.statistic = "F")#p = 0.015

#PAIRWISE COMPARISON TESTS FOR EACH FUNCTIONAL GROUP
coral.mid <- as.data.frame(post.mid[[2]]) %>% filter(stringr::str_count(contrast, 'Total Coral') == 2)
CCA.mid <- as.data.frame(post.mid[[2]]) %>% filter(stringr::str_count(contrast, 'CCA') == 2)
UPMA.mid <- as.data.frame(post.mid[[2]]) %>% filter(stringr::str_count(contrast, 'UPMA') == 2)
TURF.mid <- as.data.frame(post.mid[[2]]) %>% filter(stringr::str_count(contrast, 'TURF') == 2)
EMA.mid <- as.data.frame(post.mid[[2]]) %>% filter(stringr::str_count(contrast, 'EMA') == 2)

#apply test correction to each functional group
coral.mid$p.adj <-  p.adjust(coral.mid$p.value, method = "BH")
CCA.mid$p.adj <-  p.adjust(CCA.mid$p.value, method = "BH")
UPMA.mid$p.adj <-  p.adjust(UPMA.mid$p.value, method = "BH") 
TURF.mid$p.adj <-  p.adjust(TURF.mid$p.value, method = "BH") 
EMA.mid$p.adj <-  p.adjust(EMA.mid$p.value, method = "BH") 

post.mid <- bind_rows(coral.mid, CCA.mid, UPMA.mid, TURF.mid, EMA.mid)
post.mid$strata <- "Mid"

#SUBSET DEEP
modR.deep<-svyglm(Percent/100 ~ New*OBS_YEAR, design=des.deep, family = quasibinomial(link = "logit"))
deep <- bind_cols("Deep", car::Anova(modR.deep, type = 3, test.statistic = "F"))
post.deep <- emmeans(modR.deep, specs = pairwise ~ OBS_YEAR*New, adjust = "none")

#Single Factor tests
#deep.coral <- car::Anova(svyglm(Percent/100 ~ OBS_YEAR, design=subset(des, DEPTH_BIN == "Deep" & New == "Total Coral"), family = quasibinomial(link = "logit")), type = 3, test.statistic = "F")#NS p = 0.11
#deep.CCA <- car::Anova(svyglm(Percent/100 ~ OBS_YEAR, design=subset(des, DEPTH_BIN == "Deep" & New == "CCA"), family = quasibinomial(link = "logit")), type = 3, test.statistic = "F")#p = 0.0005
#deep.UPMA <- car::Anova(svyglm(Percent/100 ~ OBS_YEAR, design=subset(des, DEPTH_BIN == "Deep" & New == "UPMA"), family = quasibinomial(link = "logit")), type = 3, test.statistic = "F")#p = 0.003
#deep.TURF <- car::Anova(svyglm(Percent/100 ~ OBS_YEAR, design=subset(des, DEPTH_BIN == "Deep" & New == "TURF"), family = quasibinomial(link = "logit")), type = 3, test.statistic = "F")#p = 0.0003
#deep.EMA <- car::Anova(svyglm(Percent/100 ~ OBS_YEAR, design=subset(des, DEPTH_BIN == "Deep" & New == "EMA"), family = quasibinomial(link = "logit")), type = 3, test.statistic = "F")#p < 0.0001

#PAIRWISE COMPARISON TESTS FOR EACH FUNCTIONAL GROUP
coral.deep <- as.data.frame(post.deep[[2]]) %>% filter(stringr::str_count(contrast, 'Total Coral') == 2)
CCA.deep <- as.data.frame(post.deep[[2]]) %>% filter(stringr::str_count(contrast, 'CCA') == 2)
UPMA.deep <- as.data.frame(post.deep[[2]]) %>% filter(stringr::str_count(contrast, 'UPMA') == 2)
TURF.deep <- as.data.frame(post.deep[[2]]) %>% filter(stringr::str_count(contrast, 'TURF') == 2)
EMA.deep <- as.data.frame(post.deep[[2]]) %>% filter(stringr::str_count(contrast, 'EMA') == 2)

#apply test correction to each functional group
coral.deep$p.adj <-  p.adjust(coral.deep$p.value, method = "BH")
CCA.deep$p.adj <-  p.adjust(CCA.deep$p.value, method = "BH") 
UPMA.deep$p.adj <-  p.adjust(UPMA.deep$p.value, method = "BH") 
TURF.deep$p.adj <-  p.adjust(TURF.deep$p.value, method = "BH") 
EMA.deep$p.adj <-  p.adjust(EMA.deep$p.value, method = "BH") 

post.deep <- bind_rows(coral.deep, CCA.deep, UPMA.deep, TURF.deep, EMA.deep)
post.deep$strata <- "Deep"

#put all post hoc tests together for export
post.hocs <- dplyr::bind_rows(post.deep, post.mid, post.shal)

#export glm  and post hoc outputs
write.csv(post.hocs, "T:/Benthic/Projects/Swains 2023 Benthic Analysis/Tables/BenthicCover_posthoc_2010_lonelypsu.csv")
write.csv(bind_rows(shal, mid, deep), "T:/Benthic/Projects/Swains 2023 Benthic Analysis/Tables/BenthicCover_anova.csv")

des<-svydesign(id=~1, strata=~ Strat_conc, weights=~sw,data= filter(site.sw.pool, New!= "Other")) 

##########ALL TAXA BENTHIC COVER PLOT#####
#create strata level estimates of percent cover from the survey package
strata.means <- svyby(~Percent, ~Strat_conc, des, svymean) #calculate mean and standard error based on weighting
strata.means <- strata.means %>% tidyr::separate(Strat_conc, sep = "_", into = c("OBS_YEAR", "DEPTH_BIN", "TAXA")) #unconcatenate strata variables
strata.means$OBS_YEAR <- as.Date(ISOdate(strata.means$OBS_YEAR,1,1)) #convert OBS_YEAR to date


#create colors for functional groups
taxalist <- c("Total Coral", "CCA","UPMA" , "TURF","EMA")
taxacolors <- c("black", "hotpink1", "palegreen4", "limegreen","orange2")

strata.means$TAXA <- factor(strata.means$TAXA, levels = taxalist) #make taxa a factor
strata.means$DEPTH_BIN <- factor(strata.means$DEPTH_BIN, levels = c("Shallow", "Mid", "Deep")) #make depth strata a factor
panels <- data.frame(DEPTH_BIN =c("Shallow", "Mid", "Deep"), LAB = c("a", "b", "c"), x = as.Date(("2009-09-01")), y = c(50, 45, 46))#create panel labels
panels$DEPTH_BIN <- factor(panels$DEPTH_BIN, levels = c("Shallow", "Mid", "Deep")) #make panel label depths a factor for plot ordering

lines <- ggplot(strata.means, aes(x = OBS_YEAR, y = Percent, color = TAXA)) + 
  geom_point(size = 4)+ 
  geom_line()+
  geom_errorbar(aes(ymin = Percent - se, ymax = Percent + se), width = 1)+
  facet_wrap(~DEPTH_BIN, scales = "free_y", ncol = 1)+
  scale_x_date(date_labels = "%Y", 
               date_breaks = "1 year", 
               limits = as.Date(c("2009-08-01", "2023-06-01")),
               expand = c(0,0))+
  geom_text(inherit.aes = F, data = panels, aes(x = x, y = Inf, label = LAB), vjust = 1.15, hjust = -.1,size = 6)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        strip.background = element_blank(),legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        strip.text = element_text(size = 16))+
  labs(x = "Year",
       y = "Percent Cover",
       color = "Functional Group")+
  scale_color_manual(values = taxacolors)+
  geom_vline(xintercept = as.Date("2016-01-01"), color = "red", linetype = 2)


ggsave("T:/Benthic/Projects/Swains 2023 Benthic Analysis/Plots/Figure2.png", lines, height = 8, width = 11)
write.csv(strata.means,"Data/Strata_level_benthic_cover.csv")
