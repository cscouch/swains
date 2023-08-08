#Code for running PERMANOVA, SIMPER, and nMDS plots on Swains coral density data
#written by Brittany

#load libraries
library(tidyverse)
library(vegan)


rm(list=ls())



####PREP MULTIVARIATE DATA MATRIX---------------
setwd("C:/github/swains/data")
adults <- read.csv("C:/Users/Brittany.Huntington/Desktop/Swains/Data/BenthicREA_sitedata_TAXONCODE.csv")

#filter down to columns we care about
dat <- adults %>%
  filter(ISLAND == "Swains") %>%
  select(SITE, TAXONCODE, AdColDen, OBS_YEAR, DEPTH_BIN, LATITUDE, LONGITUDE, MIN_DEPTH_M, MAX_DEPTH_M)%>%
  mutate_if(is.character,as.factor) 

levels(dat$TAXONCODE) 

#drop taxa with 0 density at Swains and build species matrix
dat <- dat %>%
  filter(AdColDen > "0")%>%
  spread(key = "TAXONCODE", value = "AdColDen")%>%
  mutate_all(~replace_na(.,0))%>%
  select(-SSSS, -UNKN)

#check out colnames for merging certain 'problem taxa' and to align with classifications used in 2023
colnames(dat)

#updating TAXA CODE values to reflect our 2023 scoring and reduce our best guesses on ID errors from 2015 & 2018
dat <- dat %>%
  mutate(
    ACSP = ACSP + AVER +ARET +AHYA,
    ASTS = ASTS + ACUR,
    COSP = CEXE + CCOL + COSP,
    FUSP = FUSP + CYSP,
    LESP = LINC + LMYC +LESP,
    MOSP = MINC + MOSP,
    PMVC = PDAM + PDAN + PMVC +POCS,
    PAVS = PAVS + PDIF,
    PSSP = PSSP + PHAI,
    POSP = POSP + PSOL,
    STYS = SPIS + STYS
  ) %>%
  select(-AVER, -ARET, -AHYA, -ACUR, -CEXE, -CCOL, -CMAY, -CYSP, -LINC, -LMYC, -MINC, -PDAM, -PDAN,-POCS,  -PDIF, -PHAI, -PSOL, -SPIS,
          -CMAY, -LEPT)

###check out # of sites by year, depth.
tbl_N<- dat %>% 
  group_by(OBS_YEAR, DEPTH_BIN) %>%
  count()%>%
  spread(DEPTH_BIN, n)

tbl_N #unbalanced design; 2018 has twice the sample size across all depths as 2015


####PERMANOVA---------------

#create dataframes for (1) coral taxa & (2) driver varaiable
groups <- select (dat, SITE:MAX_DEPTH_M)
taxa <- select(dat, ACSP:STYS)

#create a heat map of taxa matrix
range(taxa) #ideally btwn 0-10
taxa <- (taxa)^(1/2) #square root tranformation to get data range between 0-10
taxa_mtrx <- as.matrix(taxa) #convert to matrix
heatmap(taxa_mtrx) #visualize distribution of dataframe
range(taxa) #ideally btwn 0-10

##OPTIONAL: eliminate rare species
taxa_pa <- decostand(taxa, "pa") #use 'decostand' function to converts matrix to presence/absence per site
taxa_sum <- apply(taxa_pa, 2, sum) #calculate sum per species over columns
sort(taxa_sum) #see which taxa are rare
taxa_reduced <- taxa[, !taxa_sum<2] # removes taxa only seen once.


#OPTIONAL: adding dummy variable column to species maxtrix if needed as can't run PERMANOVA on 'empty' sites
#taxa_reduced$dummy <- min(apply(taxa, 1, function(x) min(x[x>0]))) #set dummy value to lowest non-zero density value in taxa

#PERMANOVA 
#completely random design (no restrictions on permutations, use of residual in error term of all test statistics)

basic <- adonis2(taxa_reduced ~ OBS_YEAR + DEPTH_BIN, data = groups, permutations = 999, method = 'bray', by = 'margin')
basic

#add in interaction term (this version of vegan has the surprising action that when applied to an A*B model, outputs proceeds to ignore the main effects!)
basic2 <- adonis2(taxa_reduced ~ DEPTH_BIN * OBS_YEAR, data = groups, permutations = 999, method = 'bray', by = 'margin')
basic2

#set blocking factor of depth bin
perm <- with(groups, how(nperm = 999, blocks = DEPTH_BIN))  #set depth bin as a blocking factor
pmv <- adonis2(taxa_reduced ~ DEPTH_BIN + OBS_YEAR, data = groups, permutations = perm, method = "bray", sqrt.dist = FALSE, by = "terms")
pmv

#create a permutation object for our uniqiue StRS study design: Split Plot Analysis (Treatment term = OBS_YEAR)
#specify that sites are to be freely permuted within blocks but that blocks are not allowed to permute
perm <- how(within = Within(type = "free"),
               plots = Plots(type = "none"),
               blocks = groups$DEPTH_BIN,
               nperm = 999,
               observed = TRUE)

check(groups, control = perm) #check possible # of permutations---> lots!

pmv_complex <- adonis2(taxa_reduced  ~ DEPTH_BIN + OBS_YEAR, data = groups, permutations = perm, method = "bray", by = "terms")
pmv_complex

#DEPTH_BIN is added first to ensure that it accounts for as much of the variation as possible.  There are 2 df for 
#block, and this accounts for all of the variation partitioned to the block terms in the whole plot analysis.  
#Since the blocks were not permuted, the P-value for this term should be ignored!
#DEPTH_BIN accounts for 23% of variation in the whole plot analysis, OBS_YEAR only 7% but this is significant


#check PERMANOVA assumption of equal dispersion among groups/levels for each factor
taxa_distmat <- vegdist((taxa_reduced), method = "bray")

#Year
bd <-  betadisper(taxa_distmat, groups$OBS_YEAR)
boxplot(bd)
anova(bd) #passes assumption
#TukeyHSD(bd, ordered = FALSE, conf.level = 0.95) #optional post-hoc to see where dispersion may differ among factor levels

#DEPTH
bd <-  betadisper(taxa_distmat, groups$DEPTH_BIN)
boxplot(bd)
anova(bd) #passes assumption


#pairwise comparisons
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

pw_YEAR <- pairwise.adonis2(taxa_reduced ~ OBS_YEAR, data = groups)
pw_YEAR #all years differ from one another
#pairwise post-hoc tests

pw_DEPTH <- pairwise.adonis2(taxa_reduced ~ DEPTH_BIN, data = groups)
pw_DEPTH #All depths differ from one another---> though we used this as a blocking factor not a main effect in the PERMANOVA model

####SIMPER------------------------------
#explore which specific taxa are differing among OBS_YEAR
sim <- simper(taxa_reduced, group = groups$OBS_YEAR, permutations = perm)
sim
summary(sim) 

#explore which specific taxa are differing among DEPTH_BIN (zonation by depth)
sim <- simper(taxa_reduced, group = groups$DEPTH_BIN)
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

sig.species.scores <- subset(species.scores, pval<=0.05) #subset data to show species significant at 0.05


#plot
g1 <- ggplot(data.scores, aes(x=NMDS1, y=NMDS2))+ #sets up the plot
  geom_point(aes(NMDS1, NMDS2, colour = factor(data.scores$OBS_YEAR)), size = 2)+ #adds site points to plot, shape determined by Landuse, colour determined by Management
  coord_fixed()+
  theme_classic()+ 
  labs(colour = "YEAR")+ # add legend label
  theme(legend.position = "right", legend.text = element_text(size = 12), legend.title = element_text(size = 12), axis.text = element_text(size = 10)) # add legend at right of plot
g1

#add in species vectors for those with p-values <0.05 from the envfit
g2 <- g1+
  geom_segment(data = sig.species.scores, aes(x = 0, xend=NMDS1, y=0, yend=NMDS2), arrow = arrow(length = unit(0.25, "cm")), colour = "grey10", lwd=0.3) + #add vector arrows of significant species
  ggrepel::geom_text_repel(data = sig.species.scores, aes(x=NMDS1, y=NMDS2, label = SPECIES), cex = 3, direction = "both", segment.size = 0.25)+ #add labels for species, use ggrepel::geom_text_repel so that labels do not overlap
  labs(title = "Ordination with species vectors")
g2


# function for 95% confidence ellipses 
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

# data for labelling the ellipse
NMDS.mean=aggregate(data.scores[,c("NMDS1", "NMDS2")], 
                    list(group = data.scores$OBS_YEAR), mean)
#add ellipses to plot
g3 <- g2+ 
  geom_path(data = df_ell, aes(x = NMDS1, y = NMDS2, group = OBS_YEAR, color = OBS_YEAR)) #this is the ellipse, seperate ones by Site. 
g3


ggplot2::ggsave ("nMDS_all_years_taxa_reduced.jpeg", width = 6, height = 5, units = 'in')

####nMDS plot: ggplot depth bin---------------------------
data.scores <- as.data.frame(scores(nMDS, "site"))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores$site <- rownames(data.scores)  # create a column of site names, from the rownames of data.scores
data.scores$DEPTH_BIN <- as.factor(groups$DEPTH_BIN)  #  add the grp variable created earlier
head(data.scores)  #look at the data

species.scores <- as.data.frame(scores(nMDS, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores
species.scores <- cbind(species.scores, pval = taxa.spp.fit$vectors$pvals)
head(species.scores)

sig.species.scores <- subset(species.scores, pval<=0.05) #subset data to show species significant at 0.05


#plot
g1 <- ggplot(data.scores, aes(x=NMDS1, y=NMDS2))+ #sets up the plot
  geom_point(aes(NMDS1, NMDS2, colour = factor(data.scores$DEPTH_BIN)), size = 2)+ #adds site points to plot, shape determined by Landuse, colour determined by Management
  coord_fixed()+
  theme_classic()+ 
  labs(colour = "YEAR")+ # add legend label
  theme(legend.position = "right", legend.text = element_text(size = 12), legend.title = element_text(size = 12), axis.text = element_text(size = 10)) # add legend at right of plot
g1


g2 <- g1+
  geom_segment(data = sig.species.scores, aes(x = 0, xend=NMDS1, y=0, yend=NMDS2), arrow = arrow(length = unit(0.25, "cm")), colour = "grey10", lwd=0.3) + #add vector arrows of significant species
  ggrepel::geom_text_repel(data = sig.species.scores, aes(x=NMDS1, y=NMDS2, label = species), cex = 3, direction = "both", segment.size = 0.25) #add labels for species, use ggrepel::geom_text_repel so that labels do not overlap
g2


# function for ellipses 
veganCovEllipse <- function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

#data for ellipse, in this case using the YEAR factor
df_ell <- data.frame() #sets up a data frame before running the function.
for(g in levels(data.scores$DEPTH_BIN)){
  df_ell <- rbind(df_ell, cbind(as.data.frame(with(data.scores [data.scores$DEPTH_BIN==g,],
                                                   veganCovEllipse(cov.wt(cbind(NMDS1,NMDS2),wt=rep(1/length(NMDS1),length(NMDS1)))$cov,center=c(mean(NMDS1),mean(NMDS2))))) 
                                ,DEPTH_BIN=g))
}

# data for labelling the ellipse
NMDS.mean=aggregate(data.scores[,c("NMDS1", "NMDS2")], 
                    list(group = data.scores$DEPTH_BIN), mean)
#add ellipses to plot
g3 <- g2+ 
  geom_path(data = df_ell, aes(x = NMDS1, y = NMDS2, group = DEPTH_BIN, color = DEPTH_BIN)) #this is the ellipse, seperate ones by Site. 
g3

setwd("C:/github/swains/colony_density")
ggplot2::ggsave ("nMDS_depth_taxa_reduced.jpeg", width = 6, height = 5, units = 'in')


