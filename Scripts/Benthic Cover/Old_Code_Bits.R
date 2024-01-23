##Random bits of swians code that we don't really need

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
  dplyr::summarise(TEMP_C = mean(TEMP_C))
swains.temps.monthly$DATE <- as.Date(paste(swains.temps.monthly$YEAR, swains.temps.monthly$MONTH, "1",sep = "-"), format = "%Y-%m-%d")

colors <- c("#AA4499" ,"#44BB99", "#FF8C00")

dhw <- c(2016-03-31, 2016-06-19)

STR_temps <- 
  
  ggplot(swains.temps.monthly[swains.temps.monthly$YEAR %in% c(2015,2016),], aes(x = DATE, y = TEMP_C, color = DEPTH_BIN))+
  geom_line()+
  annotate("rect", xmin = as.Date("2016-03-31", format = "%Y-%m-%d"), xmax = as.Date("2016-06-19", format = "%Y-%m-%d"),ymin = -Inf, ymax = Inf,
           alpha = .2)+
  scale_color_manual(values = colors)+
  theme_classic()+
  labs(x = "Date",
       y = expression("Temperature " ( degree~C)),
       color = "Depth Bin")+
  theme(legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        axis.text = element_text(size = 14, color = "black"),
        axis.title = element_text(size = 16),
        axis.title.x = element_blank())

ggsave("T:/Benthic/Projects/Swains 2023 Benthic Analysis/plots/STR_TEMPS.png", STR_temps, height = 8, width = 8)


######Calcifiers#####
site.sw.calc <- site.sw
site.sw.calc$New <- recode_factor(site.sw.calc$New, "Total Coral" = "CALC", .default = "NC")
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