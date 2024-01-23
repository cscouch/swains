#Code for running PERMANOVA, SIMPER, and nMDS plots on Swains coral density data
#written by Brittany and modified by Mia for Swains project

#load libraries
library(tidyverse)
library(readr)
library(vegan)
library(dplyr)

rm(list=ls())
dir = Sys.info()[7]
setwd(paste0("C:/Users/", dir, "/Documents/GitHub/swains/"))

####PREP MULTIVARIATE DATA MATRIX---------------
adults <- read.csv("Data/Swains_sitedata_TAXONCODE_MORPH.csv")

#filter down to columns we care about
dat <- adults %>%
  filter(ISLAND == "Swains") %>%
  dplyr::select(SITE, TAXONCODE_2, AdColDen, OBS_YEAR, DEPTH_BIN)%>%
  mutate_if(is.character,as.factor) 

levels(dat$TAXONCODE_2) 

#drop taxa with 0 density at Swains, filter by depth bins, and build species matrix per bin
dat_shallow <- dat %>%
  filter(AdColDen > "0")%>%
  filter(DEPTH_BIN == "Shallow") %>%
  spread(key = "TAXONCODE_2", value = "AdColDen")%>%
  mutate_all(~replace_na(.,0))

dat_mid <- dat %>%
  filter(AdColDen > "0")%>%
  filter(DEPTH_BIN == "Mid") %>%
  spread(key = "TAXONCODE_2", value = "AdColDen")%>%
  mutate_all(~replace_na(.,0))

dat_deep <- dat %>%
  filter(AdColDen > "0")%>%
  filter(DEPTH_BIN == "Deep") %>%
  spread(key = "TAXONCODE_2", value = "AdColDen")%>%
  mutate_all(~replace_na(.,0))

###check out # of sites by year, depth.
tbl_N<- dat %>% 
  group_by(OBS_YEAR, DEPTH_BIN) %>%
  count()%>%
  spread(DEPTH_BIN, n)
tbl_N #unbalanced design; 2018 has twice the sample size across all depths as 2015


##########################################################################################################
####PERMANOVA and nMDS for Shallow ----------------------------------------------------------------------

#create dataframes for (1) coral taxa & (2) driver varaiable
groups <- select (dat_shallow, SITE:DEPTH_BIN)
taxa <- select(dat_shallow, ACSP_BR:SPIS)

#create a heat map of taxa matrix
range(taxa) #ideally btwn 0-10
taxa <- (taxa)^(1/2) #square root tranformation to get data range between 0-10
taxa_mtrx <- as.matrix(taxa) #convert to matrix
heatmap(taxa_mtrx) #visualize distribution of dataframe
range(taxa) #ideally btwn 0-10

#eliminate rare species
taxa_pa <- decostand(taxa, "pa") #use 'decostand' function to converts matrix to presence/absence per site
taxa_sum <- apply(taxa_pa, 2, sum) #calculate sum per species over columns
sort(taxa_sum) #view taxa prevalence
taxa_reduced <- taxa[, !taxa_sum<2] # removes taxa only seen once.

#OPTIONAL: adding dummy variable column to species maxtrix if needed as can't run PERMANOVA on 'empty' sites
#taxa_reduced$dummy <- min(apply(taxa, 1, function(x) min(x[x>0]))) #set dummy value to lowest non-zero density value in taxa

#PERMANOVA 
pmv_shallow <- adonis2(taxa_reduced  ~  OBS_YEAR, data = groups, permutations = 999, method = "bray", by = "terms")
pmv_shallow 

#pairwise post-hoc tests to see which years differed from one another
#load function
pairwise.adonis2 <- function(x, data, strata = NULL, nperm=999, ... ) {
  
  #describe parent call function
  ststri <- ifelse(is.null(strata),'Null',strata)
  fostri <- as.character(x)
  #list to store results
  
  #copy model formula
  x1 <- x
  # extract left hand side of formula
  lhs <- eval(x1[[2]], environment(x1), globalenv())
  environment(x1) <- environment()
  # extract factors on right hand side of formula
  rhs <- x1[[3]]
  # create model.frame matrix
  x1[[2]] <- NULL
  rhs.frame <- model.frame(x1, data, drop.unused.levels = TRUE)
  
  # create unique pairwise combination of factors
  co <- combn(unique(as.character(rhs.frame[,1])),2)
  
  # create names vector
  nameres <- c('parent_call')
  for (elem in 1:ncol(co)){
    nameres <- c(nameres,paste(co[1,elem],co[2,elem],sep='_vs_'))
  }
  #create results list
  res <- vector(mode="list", length=length(nameres))
  names(res) <- nameres
  
  #add parent call to res
  res['parent_call'] <- list(paste(fostri[2],fostri[1],fostri[3],', strata =',ststri, ', permutations',nperm ))
  
  
  #start iteration trough pairwise combination of factors
  for(elem in 1:ncol(co)){
    
    #reduce model elements
    if(inherits(eval(lhs),'dist')){
      xred <- as.dist(as.matrix(eval(lhs))[rhs.frame[,1] %in% c(co[1,elem],co[2,elem]),
                                           rhs.frame[,1] %in% c(co[1,elem],co[2,elem])])
    }else{
      xred <- eval(lhs)[rhs.frame[,1] %in% c(co[1,elem],co[2,elem]),]
    }
    
    mdat1 <-  data[rhs.frame[,1] %in% c(co[1,elem],co[2,elem]),]
    
    # redefine formula
    if(length(rhs) == 1){
      xnew <- as.formula(paste('xred',as.character(rhs),sep='~'))
    }else{
      xnew <- as.formula(paste('xred' ,
                               paste(rhs[-1],collapse= as.character(rhs[1])),
                               sep='~'))}
    
    #pass new formula to adonis
    if(is.null(strata)){
      ad <- adonis2(xnew,data=mdat1, ... )
    }else{
      perm <- how(nperm = nperm)
      setBlocks(perm) <- with(mdat1, mdat1[,ststri])
      ad <- adonis2(xnew,data=mdat1,permutations = perm, ... )}
    
    res[nameres[elem+1]] <- list(ad[1:5])
  }
  #names(res) <- names
  class(res) <- c("pwadstrata", "list")
  return(res)
}

pw_YEAR <- pairwise.adonis2(taxa_reduced ~ OBS_YEAR, data = groups, p.adjust="BH")
pw_YEAR #all years differ from one another

#check PERMANOVA assumption of equal dispersion among groups/levels for each factor
taxa_distmat <- vegdist((taxa_reduced), method = "bray")

#factor: Year
bd <-  betadisper(taxa_distmat, groups$OBS_YEAR)
boxplot(bd)
anova(bd) #passes assumption
#TukeyHSD(bd, ordered = FALSE, conf.level = 0.95) #optional post-hoc to see where dispersion may differ among factor levels


####SIMPER: shallow------------------------------
#explore which specific taxa are differing among OBS_YEAR
sim <- simper(taxa_reduced, group = groups$OBS_YEAR, permutations = 999)
sim
summary(sim) 


####nMDS--------------------------------
nMDS <- metaMDS(taxa_reduced, distance = "bray", k = 3, maxit = 999, trymax = 250, wascores = TRUE)
nMDS
stressplot(nMDS)

taxa.spp.fit <- envfit(nMDS, taxa_reduced, permutations = 999) # this fits species vectors


####nMDS plot: ggplot year---------------------------
data.scores <- as.data.frame(scores(nMDS, "site"))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores$SITE <- rownames(data.scores)  # create a column of site names, from the rownames of data.scores
data.scores$OBS_YEAR <- as.factor(groups$OBS_YEAR)  #  add the grp variable created earlier
head(data.scores)  #look at the data

species.scores <- as.data.frame(scores(nMDS, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores$SPECIES <- rownames(species.scores)  # create a column of species, from the rownames of species.scores
species.scores <- cbind(species.scores, pval = taxa.spp.fit$vectors$pvals)
head(species.scores)

species.scores$SPECIES <-  as.factor(species.scores$SPECIES)
#sig.species.scores <- subset(species.scores, pval<=0.05) #subset data to show species significant at 0.05
sig.species.scores <- subset(species.scores, SPECIES %in% c("MOSP_FO", "MOSP_EM", "PGWC", "PMVC", "POSP_EM")) #subset to specific taxa contributing to 50% cummulative difference

# function for ellipses 
veganCovEllipse <- function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

#data for ellipse, in this case using the YEAR factor
df_ell <- data.frame() #sets up a data frame before running the function.
for(g in levels(data.scores$OBS_YEAR)){
  df_ell <- rbind(df_ell, cbind(as.data.frame(with(data.scores [data.scores$OBS_YEAR==g,],
                                                   veganCovEllipse(cov.wt(cbind(NMDS1,NMDS2),wt=rep(1/length(NMDS1),length(NMDS1)))$cov,center=c(mean(NMDS1),mean(NMDS2))))) 
                                ,OBS_YEAR=g))
}

##########################################################################################################
####PERMANOVA and nMDS for Mid ----------------------------------------------------------------------

#create dataframes for (1) coral taxa & (2) driver varaiable
groups_m <- select (dat_mid, SITE:DEPTH_BIN)
taxa_m <- select(dat_mid, ACSP_BR:SPIS)

#create a heat map of taxa matrix
taxa_m <- (taxa_m)^(1/2) #square root tranformation to get data range between 0-10
range(taxa_m) #ideally btwn 0-10

##OPTIONAL: eliminate rare species
taxa_pa_m <- decostand(taxa_m, "pa") #use 'decostand' function to converts matrix to presence/absence per site
taxa_mum_m <- apply(taxa_pa_m, 2, sum) #calculate sum per species over columns
sort(taxa_mum_m) #see which taxa are rare
taxa_reduced_m <- taxa_m[, !taxa_mum_m<2] # removes taxa only seen once.

#PERMANOVA 
pmv_mid <- adonis2(taxa_reduced_m  ~  OBS_YEAR, data = groups_m, permutations = 999, method = "bray", by = "terms")
pmv_mid

#check PERMANOVA assumption of equal dispersion among groups/levels for each factor
taxa_distmat_m <- vegdist((taxa_reduced_m), method = "bray")

#factor: Year
bd_mid <-  betadisper(taxa_distmat_m, groups_m$OBS_YEAR)
boxplot(bd_mid)
anova(bd_mid) #passes assumption
#TukeyHSD(bd, ordered = FALSE, conf.level = 0.95) #optional post-hoc to see where dispersion may differ among factor levels


####SIMPER: mid------------------------------
#explore which specific taxa are differing among OBS_YEAR
sim <- simper(taxa_reduced_m, group = groups_m$OBS_YEAR, permutations = 999)
sim
summary(sim) 


####nMDS--------------------------------
nMDS_m <- metaMDS(taxa_reduced_m, distance = "bray", k = 3, maxit = 999, trymax = 250, wascores = TRUE)
nMDS_m
stressplot(nMDS_m)

taxa.spp.fit_m <- envfit(nMDS_m, taxa_reduced_m, permutations = 999) # this fits species vectors


####nMDS plot: ggplot year---------------------------
data.scores_m <- as.data.frame(scores(nMDS_m, "site"))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores_m$SITE <- rownames(data.scores_m)  # create a column of site names, from the rownames of data.scores
data.scores_m$OBS_YEAR <- as.factor(groups_m$OBS_YEAR)  #  add the grp variable created earlier
head(data.scores_m)  #look at the data

species.scores_m <- as.data.frame(scores(nMDS_m, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores_m$SPECIES <- rownames(species.scores_m)  # create a column of species, from the rownames of species.scores
species.scores_m <- cbind(species.scores_m, pval = taxa.spp.fit_m$vectors$pvals)
head(species.scores_m)

species.scores_m$SPECIES <-  as.factor(species.scores_m$SPECIES)
#sig.species.scores <- subset(species.scores, pval<=0.05) #subset data to show species significant at 0.05
sig.species.scores_m <- subset(species.scores_m, SPECIES %in% c("MOSP_FO", "MOSP_EM", "PMAL", "PMVC", "PSSP_EM","LMYC")) #subset to specific taxa contributing to 50,% cummulative difference

#data for ellipse, in this case using the YEAR factor
df_ell_m <- data.frame() #sets up a data frame before running the function.
for(g in levels(data.scores_m$OBS_YEAR)){
  df_ell_m <- rbind(df_ell_m, cbind(as.data.frame(with(data.scores_m [data.scores_m$OBS_YEAR==g,],
                                                       veganCovEllipse(cov.wt(cbind(NMDS1,NMDS2),wt=rep(1/length(NMDS1),length(NMDS1)))$cov,center=c(mean(NMDS1),mean(NMDS2))))) 
                                    ,OBS_YEAR=g))
}


##########################################################################################################
####PERMANOVA and nMDS for Deep ----------------------------------------------------------------------

#create dataframes for (1) coral taxa & (2) driver varaiable
groups_d <- select (dat_deep, SITE:DEPTH_BIN)
taxa_d <- select(dat_deep, ASTS_MD:SPIS)

#create a heat map of taxa matrix
taxa_d <- (taxa_d)^(1/2) #square root tranformation to get data range between 0-10
range(taxa_d) #ideally btwn 0-10

##OPTIONAL: eliminate rare species
taxa_pa_d <- decostand(taxa_d, "pa") #use 'decostand' function to converts matrix to presence/absence per site
taxa_dum_d <- apply(taxa_pa_d, 2, sum) #calculate sum per species over columns
sort(taxa_dum_d) #see which taxa are rare
taxa_reduced_d <- taxa_d[, !taxa_dum_d<2] # removes taxa only seen once.

#PERMANOVA 
pmv_deep <- adonis2(taxa_reduced_d  ~  OBS_YEAR, data = groups_d, permutations = 999, method = "bray", by = "terms")
pmv_deep

#check PERMANOVA assumption of equal dispersion among groups/levels for each factor
taxa_distmat_d <- vegdist((taxa_reduced_d), method = "bray")

#factor: Year
bd_deep <-  betadisper(taxa_distmat_d, groups_d$OBS_YEAR)
boxplot(bd_deep)
anova(bd_deep) #passes assumption
#TukeyHSD(bd, ordered = FALSE, conf.level = 0.95) #optional post-hoc to see where dispersion may differ among factor levels


####SIMPER: deep------------------------------
#explore which specific taxa are differing among OBS_YEAR
sim <- simper(taxa_reduced, group = groups$OBS_YEAR, permutations = 999)
sim
summary(sim) 


####nMDS--------------------------------
nMDS_d <- metaMDS(taxa_reduced_d, distance = "bray", k = 3, maxit = 999, trymax = 250, wascores = TRUE)
nMDS_d
stressplot(nMDS_d)

taxa.spp.fit_d <- envfit(nMDS_d, taxa_reduced_d, permutations = 999) # this fits species vectors


####nMDS plot: ggplot year---------------------------
data.scores_d <- as.data.frame(scores(nMDS_d, "site"))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores_d$SITE <- rownames(data.scores_d)  # create a column of site names, from the rownames of data.scores
data.scores_d$OBS_YEAR <- as.factor(groups_d$OBS_YEAR)  #  add the grp variable created earlier
head(data.scores_d)  #look at the data

species.scores_d <- as.data.frame(scores(nMDS_d, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores_d$SPECIES <- rownames(species.scores_d)  # create a column of species, from the rownames of species.scores
species.scores_d <- cbind(species.scores_d, pval = taxa.spp.fit_d$vectors$pvals)
head(species.scores_d)

species.scores_d$SPECIES <-  as.factor(species.scores_d$SPECIES)
#sig.species.scores <- subset(species.scores, pval<=0.05) #subset data to show species significant at 0.05
sig.species.scores_d <- subset(species.scores_d, SPECIES %in% c("MOSP_FO", "MOSP_EM", "PMVC", "PGWC","POSP_EM")) #subset to specific taxa contributing to 50,% cummulative difference

#data for ellipse, in this case using the YEAR factor
df_ell_d <- data.frame() #sets up a data frame before running the function.
for(g in levels(data.scores_d$OBS_YEAR)){
  df_ell_d <- rbind(df_ell_d, cbind(as.data.frame(with(data.scores_d [data.scores_d$OBS_YEAR==g,],
                                                       veganCovEllipse(cov.wt(cbind(NMDS1,NMDS2),wt=rep(1/length(NMDS1),length(NMDS1)))$cov,center=c(mean(NMDS1),mean(NMDS2))))) 
                                    ,OBS_YEAR=g))
}

#################################################################
#####Combine nMDs plots into one panel plot for paper----------------------------------------------------

#combine all data scores, species scores, and ellipse data for ease of making plot
data.scores$DEPTH_BIN = "Shallow"
data.scores_m$DEPTH_BIN = "Mid"
data.scores_d$DEPTH_BIN = "Deep"
species.scores$DEPTH_BIN = "Shallow"
species.scores_m$DEPTH_BIN = "Mid"
species.scores_d$DEPTH_BIN = "Deep"
df_ell$DEPTH_BIN = "Shallow"
df_ell_m$DEPTH_BIN = "Mid"
df_ell_d$DEPTH_BIN = "Deep"
sig.species.scores$DEPTH_BIN = "Shallow"
sig.species.scores_m$DEPTH_BIN = "Mid"
sig.species.scores_d$DEPTH_BIN = "Deep"

data.scores_all <- rbind(data.scores,data.scores_m,data.scores_d)%>%
  mutate_if(is.character,as.factor) 

species.scores_all <-rbind(species.scores,species.scores_m, species.scores_d) %>%
  mutate_if(is.character,as.factor) 

df_ell_all <-rbind(df_ell, df_ell_m,df_ell_d)%>%
  mutate_if(is.character,as.factor)

sig.species.scores_all <- rbind(sig.species.scores, sig.species.scores_m, sig.species.scores_d)%>%
  mutate_if(is.character,as.factor) 

sig.species.scores_all$SPECIES <- droplevels(sig.species.scores_all$SPECIES)
sig.species.scores_all$SPECIES <- recode_factor(sig.species.scores_all$SPECIES, "MOSP_EM" = "Enc. MOSP", 
                                                                                "MOSP_FO" = "Fol. MOSP",
                                                                                "POSP_MD" ="Mound. POSP",
                                                                                "POSP_EM" ="Enc. POSP",
                                                                                "PSSP_EM" = "Enc. PSSP")



#####Figure 3: make panel plot with a nMDS plot for each depth bin--------

data.scores_all$DEPTH_BIN <- factor(data.scores_all$DEPTH_BIN, levels = c("Shallow","Mid","Deep"))
g <- ggplot(data.scores_all, aes(x=NMDS1, y=NMDS2))+ #sets up the plot
  geom_point(aes(NMDS1, NMDS2, colour = factor(data.scores_all$OBS_YEAR)), size = 2)+ #adds site points to plot, shape determined by Landuse, colour determined by Management
  facet_wrap(~DEPTH_BIN,scales = "free",)+
  scale_colour_manual(values= alpha(c("#FF8C00","#44BB99","#AA4499")))+
  #coord_fixed()+
  theme_bw()+ 
  labs(colour = "YEAR")+ # add legend label
  theme(legend.position = "bottom", legend.text = element_text(size = 16), legend.title = element_text(size = 16), axis.text = element_text(size = 16), axis.title = element_text(size = 16)) +
  geom_segment(data = sig.species.scores_all, aes(x = 0, xend=NMDS1, y=0, yend=NMDS2), arrow = arrow(length = unit(0.25, "cm")), colour = "grey60", lwd=0.8) + #add vector arrows of significant species
  geom_path(data = df_ell_all, aes(x = NMDS1, y = NMDS2, group = OBS_YEAR, color = OBS_YEAR), lwd = 1.2)+
  ggrepel::geom_text_repel(data = sig.species.scores_all, aes(x=NMDS1, y=NMDS2, label = SPECIES), cex = 5, direction = "both", segment.size = 0.25) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(strip.text = element_text(size = 22))
g

ggsave(plot=g,file= "C:/Users/", dir, "/Documents/GitHub/swains/Plots/Figure3.jpg",width=14,height=7)


