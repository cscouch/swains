
# Using R version 4.1.0 (2021-05-18)

rm(list=ls())
dir = Sys.info()[7]
setwd(paste0("C:/Users/", dir, "/Documents/GitHub/swains/"))



library(tidyselect)
library(ggpubr)
library(survey)
library(emmeans)
library(dplyr)

#LOAD site-level data
swa <- read.csv("Data/Swains_sitedata_TAXONCODE_MORPH.csv")

#read in sector-area file
sectors<-read.csv("Data/Sectors-Strata-Areas.csv", stringsAsFactors=FALSE)
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

nrow(site.sw)

#Create contactenated Strata variable 
site.sw$Strat_conc<-paste(site.sw$OBS_YEAR, site.sw$DEPTH_BIN,sep = "_")
site.sw$OBS_YEAR<-as.factor(site.sw$OBS_YEAR)





#Subset data

# Use demography data only from the shallow sector
site.sw$OBS_YEAR<-as.factor(site.sw$OBS_YEAR)

site.sw.mospfo <- site.sw %>% filter(TAXONCODE_2 == "MOSP_FO") 
site.sw.mospem <- site.sw %>% filter(TAXONCODE_2 == "MOSP_EM")  
site.sw.pmvc <- site.sw %>% filter(TAXONCODE_2 == "PMVC") 
site.sw.pgwc <- site.sw %>% filter(TAXONCODE_2 == "PGWC") 
site.sw.pospem <- site.sw %>% filter(TAXONCODE_2 == "POSP_EM") 
site.sw.pospmd <- site.sw %>% filter(TAXONCODE_2 == "POSP_MD")


#Establish survey design with survey weights
des.pgwc<-svydesign(id=~1, strata=~ Strat_conc, weights=~sw,data=site.sw.pgwc)
des.pmvc<-svydesign(id=~1, strata=~ Strat_conc, weights=~sw,data=site.sw.pmvc)
des.mospfo<-svydesign(id=~1, strata=~ Strat_conc, weights=~sw,data=site.sw.mospfo)
des.mospem<-svydesign(id=~1, strata=~ Strat_conc, weights=~sw,data=site.sw.mospem)
des.pospem<-svydesign(id=~1, strata=~ Strat_conc, weights=~sw,data=site.sw.pospem)
des.pospmd<-svydesign(id=~1, strata=~ Strat_conc, weights=~sw,data=site.sw.pospmd)


#Shallow
site.sw.pmvc.sh <- site.sw.pmvc %>% filter(DEPTH_BIN == "Shallow") #
des.pmvc.sh<-svydesign(id=~1, strata=~ Strat_conc, weights=~sw,data=site.sw.pmvc.sh)
modR.pmvc <- svyglm(AdColCount ~ OBS_YEAR,  offset = TRANSECTAREA_ad, design = des.pmvc.sh, family = "quasipoisson")
car::Anova(modR.pmvc, type=3, test.statistic = "F") # p = 0.003

site.sw.pgwc.sh <- site.sw.pgwc %>% filter(DEPTH_BIN == "Shallow") #
des.pgwc.sh<-svydesign(id=~1, strata=~ Strat_conc, weights=~sw,data=site.sw.pgwc.sh)
modR.pgwc <- svyglm(AdColCount ~ OBS_YEAR,  offset = TRANSECTAREA_ad, design = des.pgwc.sh, family = "quasipoisson")
car::Anova(modR.pgwc, type=3, test.statistic = "F") # p < 0.00001

site.sw.mospem.sh <- site.sw.mospem %>% filter(DEPTH_BIN == "Shallow") #
des.mospem.sh<-svydesign(id=~1, strata=~ Strat_conc, weights=~sw,data=site.sw.mospem.sh)
modR.mospem <- svyglm(AdColCount ~ OBS_YEAR,  offset = TRANSECTAREA_ad, design = des.mospem.sh, family = "quasipoisson")
car::Anova(modR.mospem, type=3, test.statistic = "F") # p = 0.05498

site.sw.mospfo.sh <- site.sw.mospfo %>% filter(DEPTH_BIN == "Shallow") #
des.mospfo.sh<-svydesign(id=~1, strata=~ Strat_conc, weights=~sw,data=site.sw.mospfo.sh)
modR.mospfo <- svyglm(AdColCount ~ OBS_YEAR,  offset = TRANSECTAREA_ad, design = des.mospfo.sh, family = "quasipoisson")
car::Anova(modR.mospfo, type=3, test.statistic = "F") #p = 0.445

site.sw.pospmd.sh <- site.sw.pospmd %>% filter(DEPTH_BIN == "Shallow") #
des.pospmd.sh<-svydesign(id=~1, strata=~ Strat_conc, weights=~sw,data=site.sw.pospmd.sh)
modR.pospmd <- svyglm(AdColCount ~ OBS_YEAR,  offset = TRANSECTAREA_ad, design = des.pospmd.sh, family = "quasipoisson")
car::Anova(modR.pospmd, type=3, test.statistic = "F") #p <0.0001

site.sw.pospem.sh <- site.sw.pospem %>% filter(DEPTH_BIN == "Shallow") #
des.pospem.sh<-svydesign(id=~1, strata=~ Strat_conc, weights=~sw,data=site.sw.pospem.sh)
modR.pospem <- svyglm(AdColCount ~ OBS_YEAR,  offset = TRANSECTAREA_ad, design = des.pospem.sh, family = "quasipoisson")
car::Anova(modR.pospem, type=3, test.statistic = "F") #p 0.1046



# Post-hoc
tot.pmvc <- emmeans(modR.pmvc, specs = pairwise~OBS_YEAR, adjust = "none") # 2018 is sig less than 2015 and 2023
tot.pmvc <- emmeans(modR.pmvc,pairwise~OBS_YEAR) 
tot.pmvc.df <- as.data.frame(tot.pmvc[[2]])



# Apply test corrections
tot.pmvc.df$p.adj <- p.adjust(tot.pmvc.df$p.value, method = "BH"); tot.pmvc.df # 2018 no longer smaller than 2023

# 
# l<-c("2015","2018","2023") # p value for non-parametric tests
# ps<-matrix(NA,3,3)
# dimnames(ps)<-list(l,l)
# for(i in 1:2){
#   for(j in (i+1):3){
#     ps[i,j]<-svyranktest(AdColDen~OBS_YEAR, subset(des.pospmd, OBS_YEAR %in% l[c(i,j)]))$p.value
#   }
# }
# ps
# 
# for(i in 1:2){
#   for(j in (i+1):3){
#     ps[i,j]<-svyranktest(AdColDen~OBS_YEAR, subset(des.pgwc, OBS_YEAR %in% l[c(i,j)]))$p.value
#   }
# }
# ps
# 
# p.adjust(ps[!is.na(ps)],method="hochberg")

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
  scale_fill_manual(values = alpha(c("#FF8C00","#44BB99","#AA4499"))) +
  geom_errorbar(data = depth_mean_tot %>% filter(DEPTH_BIN == "Shallow"), 
                      aes(ymin = AdColDen-se, ymax = AdColDen+se),
                      width = .3) +
        geom_text(aes(x=OBS_YEAR,y=AdColDen+se,label=sig, group = Genus),
                  position = position_dodge(), size = 4,
                  vjust = -0.5) +
        facet_wrap(~Genus,nrow = 2, labeller = label_wrap_gen(width = 16)) +
        guides(fill="none") +
        theme_bw() +
        ylim(0,8.9) +
        theme(panel.grid.minor = element_blank(),
              panel.grid.major = element_blank())+
        ylab("Mean Adult Colonies m - 2") +
        xlab("Year")
        
        # levels = c("P. meandrina/verrucosa complex","Mounding Porites spp.", "Foliose Montipora spp.",
        #            "P. grandis/woodjonsei complex","Encrusting Porites spp.", "Encrusting Montipora spp.")),


ggsave(plot=s,file="/Plots/Figure5.jpg",width=8,height=5)

