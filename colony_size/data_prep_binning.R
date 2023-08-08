#Binning three main taxa for generating univariate response variables of size frequency distributions


#load libraries
library(tidyverse)
library(hrbrthemes)
library(viridis)
library(ggridges)


rm(list=ls())
####PREP SIZE STRUCTURE DATA---------------
setwd("C:/github/swains/data")
adults <- read.csv("CoralBelt_Adults_raw_CLEANED.csv")%>% mutate_if(is.character,as.factor) %>% filter(ISLANDCODE == "SWA") %>% droplevels()
juvs <- read.csv("CoralBelt_Juveniles_raw_CLEANED_v2.csv") %>%mutate_if(is.character,as.factor) %>% filter (ISLANDCODE == "SWA") %>% droplevels()

#filtering
adults <- adults %>% filter (TRANSECT != "2")%>%  #remove transect 2 from adults
  mutate(Fragment = ifelse(is.na(Fragment), 0, Fragment))%>%
  filter(Fragment != "-1") #Fix how fragments are recorded across years and exclude
juvs <- juvs %>% filter (TRANSECT != "2") #retaining transect 1 (and 3 and 4 from 2015)
  

adults2 <- adults[, c(5:8,12,14:17,22,42,44,47, 52)] 
juvs2 <- juvs [, c(5:9, 14:17,21,22, 24, 27,32)]

colnames(adults2)
colnames(juvs2)

adults2$OBS_YEAR <- as.factor(adults2$OBS_YEAR)
juvs2$OBS_YEAR <- as.factor(juvs2$OBS_YEAR)

######Mean and SD for adult taxa per site-------------------
dat.mean <- adults2 %>%
  filter(GENUS_CODE %in% c("MOSP", "POCS", "POSP")) %>% 
  mutate(log10.CL = log10(COLONYLENGTH))%>%
  group_by(OBS_YEAR, DEPTH_BIN, SITE, GENUS_CODE)%>%
  summarize(MEAN=mean(log10.CL)) %>%
  spread(GENUS_CODE, MEAN)
  

dat.sd <- adults2 %>%
  filter(GENUS_CODE %in% c("MOSP", "POCS", "POSP")) %>% 
  mutate(log10.CL = log10(COLONYLENGTH))%>%
  group_by(OBS_YEAR, DEPTH_BIN, SITE, GENUS_CODE)%>%
  summarize(SD=sd(log10.CL)) %>%
  spread(GENUS_CODE, SD)

dat.summary  <- left_join(dat.mean, dat.sd, by = c("OBS_YEAR", "DEPTH_BIN", "SITE"))%>%
  rename_with(~ sub(".x", ".MEAN", .x), everything())%>%
  rename_with(~ sub(".y", ".SD", .x), everything())

setwd("C:/github/swains/colony_size")
write.csv(dat.summary, "Adult_mean_SD.csv", row.names = FALSE)

######BINNING---------------
#Set to taxa
genus <- "POCS"
#genus <- "MOSP"
#genus <- "POSP"

adults2 <- adults2 %>%
  filter(GENUS_CODE %in% c(genus)) %>% #adjust this to taxa of interest
  mutate(log10.CL = log10(COLONYLENGTH))%>%
  droplevels()

#identify quantiles you want as cut offs and then cut and bin each colony into those quantiles
y <- data.frame(quantile(adults2$log10.CL, probs =  c(25,50,75)/100, na.rm = FALSE, names = TRUE, type = 9, digits = 4))
adults2$QUARTILES=cut(adults2$log10.CL,c(0,y[1,1],y[2,1],y[3,1],Inf),labels=c('Q25','Q50','Q75', "Q100"))


z <- data.frame(quantile(adults2$log10.CL, probs =  c(10,90)/100, na.rm = FALSE, names = TRUE, type = 9, digits = 4))
adults2$TAIL_BINS=cut(adults2$log10.CL,c(0,z[1,1],z[2,1],Inf),labels=c('Q10','QMed','Q90'))


#set juvenile data as bin and merge with adult data
col_order <-  c("ISLANDCODE","REEF_ZONE", "DEPTH_BIN","OBS_YEAR","SITE","TRANSECT","SEGMENT","SEGWIDTH","SEGLENGTH","COLONYLENGTH", "GENUS_CODE","TAXONNAME", "SEC_NAME",  
"TAXONCODE") #set column order to align with adults
juvs2 <- juvs2[, col_order]
juvs2 <- juvs2 %>%
  filter(GENUS_CODE %in% c(genus)) %>% #adjust this to taxa of interest
  mutate(log10.CL = log10(COLONYLENGTH), QUARTILES = "QJuv", TAIL_BINS = "QJuv")%>%
  droplevels()

dat <- rbind(adults2, juvs2)

n_distinct(dat$SITE) #check how many sites recorded this taxa

#build density in each quantile per site
site <- dat %>%
  group_by(DEPTH_BIN, OBS_YEAR, REEF_ZONE, SEC_NAME, SITE,SEGMENT)%>%
  mutate(SEG.AREA = SEGWIDTH* SEGLENGTH)%>%
  summarise(S.AREA = unique(SEG.AREA))%>%
  group_by(DEPTH_BIN, OBS_YEAR, REEF_ZONE, SEC_NAME, SITE)%>%
  summarise(SITE.AREA = sum(S.AREA))
  
quartiles <- count(dat, SITE, QUARTILES) %>%
  spread(QUARTILES, n)%>%
  mutate_all(~replace_na(.,0))

tailbins <- count(dat, SITE, TAIL_BINS) %>%
  spread(TAIL_BINS, n)%>%
  mutate_all(~replace_na(.,0))

taxa <- left_join(site, quartiles, by = "SITE", )%>%
  mutate(QJuv.D = QJuv/SITE.AREA,
         Q25.D = Q25/SITE.AREA, 
         Q50.D = Q50/SITE.AREA,
         Q75.D = Q75/SITE.AREA,
         Q100.D = Q100/SITE.AREA,
         TOTAL.D = QJuv.D + Q25.D+ Q50.D + Q75.D + Q100.D,
         QJuv.R = QJuv.D / TOTAL.D,
         Q25.R = Q25.D / TOTAL.D,
         Q50.R = Q50.D / TOTAL.D,
         Q75.R = Q75.D / TOTAL.D,
         Q100.R = Q100.D / TOTAL.D)

taxa2 <- left_join(site, tailbins, by = "SITE")%>%
  mutate(QJuv.D = QJuv/SITE.AREA,
         Q10.D = Q10/SITE.AREA, 
         QMed.D = QMed/SITE.AREA,
         Q90.D = Q90/SITE.AREA,
         TOTAL.D = QJuv.D + Q10.D+ QMed.D + Q90.D,
         QJuv.R = QJuv.D / TOTAL.D,
         Q10.R = Q10.D / TOTAL.D,
         QMed.R = QMed.D / TOTAL.D,
         Q90.R = Q90.D / TOTAL.D)


####write out data ---------
setwd("C:/github/swains/colony_size")
#write.csv(taxa, "MOSP_binned_quantiles.csv") #name to appropriate taxa
#write.csv(taxa2, "MOSP_binned_tail_bins.csv") #name to appropriate taxa


#### VISUALIZE DISTRIBUTIONS--------------------------


##plot quantile binned data
summary <- taxa %>% 
  select(DEPTH_BIN, OBS_YEAR, QJuv.R:Q100.R)%>%
  gather(key= "BIN", value = "PROP", QJuv.R, Q25.R, Q50.R, Q75.R, Q100.R )%>%
  group_by (DEPTH_BIN, OBS_YEAR, BIN) %>%
  summarise_each (funs(mean(., na.rm=T), se = sd(., na.rm=T)/sqrt(sum(!is.na(.))),), PROP)

summary$BIN <- as.factor(summary$BIN)
summary$OBS_YEAR <- as.factor(summary$OBS_YEAR)
summary$BIN <-ordered(summary$BIN, level = c("QJuv.R", "Q25.R", "Q50.R", "Q75.R", "Q100.R"))


#bar plot of mean size
g2 <- ggplot(summary, aes(x = BIN, y = mean, fill = OBS_YEAR)) + 
  geom_bar(position="dodge", stat="identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2,position=position_dodge(.9)) +
   facet_wrap("DEPTH_BIN",  ncol =1, scales = "fixed")+
  scale_fill_viridis_d() +
  theme_ipsum() +
  ylab("Proportion")+
  xlab("Size Bin")+
  ggtitle("MOSP") #name to appropriate taxa
g2

jpeg(width = 450, height = 650, filename = "MOSP quantile binned facet.jpeg") #adjust name to genera
g2
dev.off()


#plot tail bined data
summary2 <- taxa2 %>% #change to taxa
  select(DEPTH_BIN, OBS_YEAR, QJuv.R:Q90.R)%>%
  gather(key= "BIN", value = "PROP", QJuv.R, Q10.R, QMed.R, Q90.R)%>%
  group_by (DEPTH_BIN, OBS_YEAR, BIN) %>%
  summarise_each (funs(mean(., na.rm=T), se = sd(., na.rm=T)/sqrt(sum(!is.na(.))),), PROP)

summary2$BIN <- as.factor(summary2$BIN)
summary2$OBS_YEAR <- as.factor(summary2$OBS_YEAR)
summary2$BIN <-ordered(summary2$BIN, level = c("QJuv.R", "Q10.R", "QMed.R", "Q90.R"))


#bar plot of mean size
g3 <- ggplot(summary2, aes(x = BIN, y = mean, fill = OBS_YEAR)) + 
  geom_bar(position="dodge", stat="identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2,position=position_dodge(.9)) +
  facet_wrap("DEPTH_BIN",  ncol =1, scales = "fixed")+
  scale_fill_viridis_d() +
  theme_ipsum() +
  ylab("Proportion")+
  xlab("Size Bin")+
  ggtitle("MOSP") #name to appropriate taxa
g3


jpeg(width = 450, height = 650, filename = "MOSP tail binned facet.jpeg") #adjust name to genera
g3
dev.off()



