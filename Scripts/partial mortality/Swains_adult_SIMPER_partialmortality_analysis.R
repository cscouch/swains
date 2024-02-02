#Data do not include FUSP

# Using R version 4.1.0 (2021-05-18)

rm(list=ls())

library(dplyr)
library(tidyr)
library(survey)
library(multcomp)
library(emmeans)


#LOAD site-level data
swa <- read.csv("Data/Swains_sitedata_TAXONCODE_MORPH.csv")



#read in sector-area file
sectors<-read.csv("Data/Sectors-Strata-Areas.csv", stringsAsFactors=FALSE)
swa_sa<-filter(sectors,ISLAND=="Swains")
swa_sa$DEPTH_BIN<-as.factor(swa_sa$DEPTH_BIN)

NH <- swa_sa %>%
  group_by(SEC_NAME, DEPTH_BIN)%>%
  summarize(unique = NH)%>%
  group_by(DEPTH_BIN)%>%
  summarise(NH = sum(unique))

swa<-left_join(swa,NH) #merge demography data with new NH values pooled across the 2 swains sectors


#Calculate survey weights (inverse proportion weighting)
w.df<-swa %>%
  group_by(OBS_YEAR,DEPTH_BIN,NH) %>%
  summarise(n = length(unique(SITE)))

w.df$sw<-w.df$NH/w.df$n #calculate survey weights for each site

site.sw<-left_join(swa,w.df) #merge weights with site-level data
head(site.sw)



#Create contactenated Strata variable 
site.sw$Strat_conc<-paste(site.sw$OBS_YEAR, site.sw$DEPTH_BIN,sep = "_")
site.sw$OBS_YEAR<-as.factor(site.sw$OBS_YEAR)


#Subset data
site.sw$OBS_YEAR<-as.factor(site.sw$OBS_YEAR)
site.sw$DEPTH_BIN<-as.factor(site.sw$DEPTH_BIN)
# site.sw[is.na(site.sw)] <- 0


#### Partial Mortality using SIMPER results ####

# Merge demo data with new NH values pooled across the 2 swains sectors
site.sw <- filter(site.sw, DEPTH_BIN == "Shallow") %>% drop_na()


#Create concatenated Strata variable 
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
with(site.sw.pmvc, tapply((Ave.od), OBS_YEAR, shapiro.test)) # passed shapiro test
bartlett.test(Ave.od ~ Strat_conc,site.sw.pmvc) 
with(site.sw.pgwc, tapply((Ave.od), OBS_YEAR, shapiro.test)) # passed (NA for 2018)
bartlett.test(Ave.od ~ Strat_conc,site.sw.pgwc) 
with(site.sw.mospem, tapply((Ave.od), OBS_YEAR, shapiro.test)) # passed both
bartlett.test(Ave.od ~ Strat_conc,site.sw.mospem)
with(site.sw.mospfo, tapply((Ave.od), OBS_YEAR, shapiro.test)) # passed both
bartlett.test(Ave.od ~ Strat_conc,site.sw.mospfo)
with(site.sw.pospem, tapply((Ave.od), OBS_YEAR, shapiro.test)) # passed bartlett test
bartlett.test(Ave.od ~ Strat_conc,site.sw.pospem)
with(site.sw.pospmd, tapply((Ave.od), OBS_YEAR, shapiro.test)) # passed (NA for 2018)
bartlett.test(Ave.od ~ Strat_conc,site.sw.pospmd)


# Non-parametric version of svyglm
svyranktest(Ave.od ~ OBS_YEAR, design=des.pmvc, test=("KruskalWallis")) # p < 0.001
svyranktest(Ave.od ~ OBS_YEAR, design=des.pgwc, test=("KruskalWallis")) # NS
svyranktest(Ave.od ~ OBS_YEAR, design=des.pospem,DEPTH_BIN=="Shallow", test=("KruskalWallis")) # NS
svyranktest(Ave.od ~ OBS_YEAR, design=des.pospmd, test=("KruskalWallis")) # NS


#  Test fixed effects of year
modR.pgwc <- svyglm(Ave.od ~ OBS_YEAR,  design = des.pgwc)
car::Anova(modR.pgwc, type=3, test.statistic = "F") # NS

modR.mospem <- svyglm(Ave.od ~ OBS_YEAR,  design = des.mospem)
car::Anova(modR.mospem, type=3, test.statistic = "F") # NS

modR.mospfo <- svyglm(Ave.od ~ OBS_YEAR, design = des.mospfo)
car::Anova(modR.mospfo, type=3, test.statistic = "F") # NS

modR.pospmd<- svyglm(Ave.od ~ OBS_YEAR,  design = des.pospmd)
car::Anova(modR.pospmd, type=3, test.statistic = "F") # NS


# Post-hoc
tot.pmvc <- emmeans(modR.pmvc, specs = pairwise~OBS_YEAR, adjust = "none") # 2015 is sig less than 2018 and 2023
tot.pmvc.df <- as.data.frame(tot.pmvc[[2]])


# Apply test corrections
tot.pmvc.df$p.adj <- p.adjust(tot.pmvc.df$p.value, method = "BH"); tot.pmvc.df # above holds true


# Non parametric post hoc and test correction
l<-c("2015","2018","2023") # p value for non-parametric tests
ps<-matrix(NA,3,3)
dimnames(ps)<-list(l,l)
for(i in 1:2){
  for(j in (i+1):3){
    ps[i,j]<-svyranktest(Ave.od~OBS_YEAR, subset(des.pmvc, OBS_YEAR %in% l[c(i,j)]))$p.value
  }
}
ps # 2018 sig different than 2015 and 2023
p.adjust(ps[!is.na(ps)],method="hochberg") # holds true




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
od_mean_pgwc <- od_mean_pgwc %>% add_row(OBS_YEAR = "2018", DEPTH_BIN = "Shallow", Ave.od = 0, se = 0)
od_mean_pgwc$Genus <- paste("P. grandis/woodjonsei complex")
od_mean_pgwc$sig <- c(NA,NA,NA)
od_mean_pospmd <- od_mean_pospmd %>% add_row(OBS_YEAR = "2018", DEPTH_BIN = "Shallow", Ave.od = 0, se = 0)
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



# Plot Partial mortality data

s <- ggplot(od_mean_tot ,  aes(x=OBS_YEAR, y=Ave.od, fill=OBS_YEAR, group = Genus)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = alpha(c("#FF8C00","#44BB99","#AA4499"))) +
  geom_errorbar(data = od_mean_tot,aes(ymin=Ave.od-se, ymax=Ave.od+se), width = .1) +
  geom_text(aes(x=OBS_YEAR,y=Ave.od+se+3,label=sig, group = Genus),
            position = position_dodge()) +
  facet_wrap(~Genus,nrow = 2, labeller = label_wrap_gen(width = 18)) +
  guides(fill="none") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        strip.placement = "outside",
        panel.spacing = unit(0, "lines")) +
  labs(x="Year",y="Average Old Dead Extent %") +
  ylim(0,39)

ggsave(plot=s,file="Plots/FigureS5.jpg",width=8,height=5)


