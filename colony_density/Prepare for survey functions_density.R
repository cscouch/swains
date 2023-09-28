# This script will reads in the CLEANED/Analysis ready data that was generated using the following script
#C:\Users\Courtney.S.Couch\Documents\GitHub\Benthic-Scripts\REA_CoralDemography\Generate REA data\REA Coral Demography_DataPrep.R
#The script does some final tweaks to the data then generates Site-level data
#These data only include surveys conducted between 2013-2019


rm(list=ls())

#Set Run Flags
DEBUG=TRUE

#LOAD LIBRARY FUNCTIONS ...
source("C:/Users/Courtney.S.Couch/Documents/GitHub/Benthic-Scripts/Functions/Benthic_Functions_newApp_vTAOfork.R")
source("C:/Users/Courtney.S.Couch/Documents/GitHub/fish-paste/lib/core_functions.R")
source("C:/Users/Courtney.S.Couch/Documents/GitHub/fish-paste/lib/GIS_functions.R")

## LOAD benthic data
awd<-read.csv("T:/Benthic/Data/REA Coral Demography & Cover/Analysis Ready Raw data/CoralBelt_Adults_raw_CLEANED.csv")
jwd<-read.csv("T:/Benthic/Data/REA Coral Demography & Cover/Analysis Ready Raw data/CoralBelt_Juveniles_raw_CLEANED.csv")

awd=subset(awd,ISLAND=="Swains")
jwd=subset(jwd,ISLAND=="Swains")

#Final Tweaks before calculating site-level data-------------------------------------------------
#Colony fragments will be removed when you generate the site level data
#We have denoted fragments differently over the course of our dataset. This code makes sure Fragment is 0 or -1 (-1 indicates it's a fragment)

awd$Fragment<-ifelse(awd$OBS_YEAR <2018 & awd$COLONYLENGTH <5 & awd$S_ORDER=="Scleractinia",-1,awd$Fragment)
head(subset(awd,Fragment==-1& OBS_YEAR<2018)) #double check that pre 2018 fragments create
awd$Fragment[is.na(awd$Fragment)] <- 0
jwd$Fragment <- 0 # you need to add this column so that you can use the site level functions correctly
awd$METHOD<-"DIVER"
jwd$METHOD<-"DIVER"

#remove FUSP
jwd<-subset(jwd,GENUS_CODE !="FUSP")
awd<-subset(awd,GENUS_CODE !="FUSP")

#Remove colonies <1cm
jwd<-subset(jwd,COLONYLENGTH >=1)


#Quantify number of colonies between 5-10cm to add to juvenile density
awd5_10<-subset(awd,COLONYLENGTH <=10)
jcols<-colnames(jwd)
awd5_10<-subset(awd5_10,select=jcols)

#Simplify Bleaching Severity categories: in 2019 the team decided to simplify the bleaching severity from 1-5 to 1-3 to improve consistency in severity values
#This code converts the severity data collected prior to 2019 to a 1-3 scale
awd$DATE_ <- ymd(awd$DATE_)
jwd$DATE_ <- ymd(jwd$DATE_)
awd5_10$DATE_ <- ymd(awd5_10$DATE_)

#Create a look a table of all of the colony attributes- you will need this the functions below
SURVEY_COL<-c("METHOD","DATE_","SITEVISITID", "OBS_YEAR", "REGION", "REGION_NAME", "ISLAND","ISLANDCODE","SEC_NAME", "SITE", "REEF_ZONE",
              "DEPTH_BIN", "LATITUDE", "LONGITUDE","MIN_DEPTH_M","MAX_DEPTH_M","TRANSECT","SEGMENT","COLONYID","GENUS_CODE","TAXONCODE","SPCODE","COLONYLENGTH")
survey_colony<-unique(awd[,SURVEY_COL])

SURVEY_SITE<-c("METHOD","MISSIONID","DATE_","SITEVISITID", "ANALYSIS_YEAR","OBS_YEAR", "REGION", "REGION_NAME", "ISLAND","ISLANDCODE","SEC_NAME", "SITE", "REEF_ZONE",
               "DEPTH_BIN", "LATITUDE", "LONGITUDE","MIN_DEPTH_M","MAX_DEPTH_M")
survey_siteAd<-unique(awd[,SURVEY_SITE])

SURVEY_SITE<-c("METHOD","MISSIONID","DATE_","SITEVISITID", "ANALYSIS_YEAR","OBS_YEAR", "REGION", "REGION_NAME", "ISLAND","ISLANDCODE","SEC_NAME", "SITE", "REEF_ZONE",
               "DEPTH_BIN", "LATITUDE", "LONGITUDE","MIN_DEPTH_M","MAX_DEPTH_M")
survey_siteJ<-unique(jwd[,SURVEY_SITE])

write.csv(survey_siteAd,"surveysite.csv")

#We did juvenile only surveys in 2017 in PRIA, this will make sure the SV table has both adult and juv sites.
survey_site<-full_join(survey_siteJ,survey_siteAd,by = c("METHOD","MISSIONID","DATE_","SITEVISITID", "ANALYSIS_YEAR","OBS_YEAR", "REGION", "REGION_NAME", "ISLAND","ISLANDCODE","SEC_NAME", "SITE", "REEF_ZONE",
                                                         "DEPTH_BIN", "LATITUDE", "LONGITUDE","MIN_DEPTH_M","MAX_DEPTH_M"));nrow(survey_site)

survey_site<-survey_site[!duplicated(survey_site[,4]),]


# GENERATE SUMMARY METRICS at the transect-leveL BY GENUS--------------------------------------------------
#Calc_ColDen_Transect
acd.gen<-Calc_ColDen_Transect(data = awd,grouping_field = "GENUS_CODE");colnames(acd.gen)[colnames(acd.gen)=="ColCount"]<-"AdColCount";colnames(acd.gen)[colnames(acd.gen)=="ColDen"]<-"AdColDen";colnames(acd.gen)[colnames(acd.gen)=="TRANSECTAREA"]<-"TRANSECTAREA_ad"# calculate density at genus level as well as total
jcd.gen<-Calc_ColDen_Transect(jwd,"GENUS_CODE"); colnames(jcd.gen)[colnames(jcd.gen)=="ColCount"]<-"JuvColCount";colnames(jcd.gen)[colnames(jcd.gen)=="ColDen"]<-"JuvColDen";colnames(jcd.gen)[colnames(jcd.gen)=="TRANSECTAREA"]<-"TRANSECTAREA_j"
awd5_10.gen<-Calc_ColDen_Transect(awd5_10,"GENUS_CODE"); colnames(awd5_10.gen)[colnames(awd5_10.gen)=="ColCount"]<-"Ad5_10ColCount";colnames(awd5_10.gen)[colnames(awd5_10.gen)=="ColDen"]<-"Ad5_10ColDen"
awd5_10.gen<-subset(awd5_10.gen,select=-c(TRANSECTAREA))

#Calc_ColMetric_Transect
cl.gen<-Calc_ColMetric_Transect(data = awd,grouping_field = "GENUS_CODE",pool_fields = "COLONYLENGTH"); colnames(cl.gen)[colnames(cl.gen)=="Ave.y"]<-"Ave.cl" #Average % old dead
od.gen<-Calc_ColMetric_Transect(data = awd,grouping_field = "GENUS_CODE",pool_fields = "OLDDEAD"); colnames(od.gen)[colnames(od.gen)=="Ave.y"]<-"Ave.od" #Average % old dead
rd.gen<-Calc_ColMetric_Transect(data = awd,grouping_field = "GENUS_CODE",pool_fields = c("RDEXTENT1", "RDEXTENT2","RDEXTENT3")); colnames(rd.gen)[colnames(rd.gen)=="Ave.y"]<-"Ave.rd" #Average % recent dead


#CHANGE TRANSECT NUMBERS FOR JUVENILES (pre-2018 we used 3 and 4)
jcd.gen$TRANSECT[jcd.gen$TRANSECT==3]<-1
jcd.gen$TRANSECT[jcd.gen$TRANSECT==4]<-2

#Remove METHOD from dataframes before merging
acd.gen<-subset(acd.gen,select=-c(METHOD))
jcd.gen<-subset(jcd.gen,select=-c(METHOD))
awd5_10.gen<-subset(awd5_10.gen,select=-c(METHOD))
cl.gen<-subset(cl.gen,select=-c(METHOD))
od.gen<-subset(od.gen,select=-c(METHOD))
rd.gen<-subset(rd.gen,select=-c(METHOD))

#Merge density and partial mortality data together.You will need to replace the DUMMY field with the one you want
MyMerge <- function(x, y){
  df <- merge(x, y, by= c("SITE","SITEVISITID","TRANSECT","GENUS_CODE"), all.x= TRUE, all.y= TRUE)
  return(df)
}
data.gen<-Reduce(MyMerge, list(acd.gen,awd5_10.gen,jcd.gen,cl.gen,od.gen,rd.gen))
head(data.gen)

#Add METHOD back in
data.gen$METHOD<-"DIVER"

head(data.gen)

#Change NaN to NA
is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))

data.gen[is.nan(data.gen)] <- NA

#There will be some NAs when you merge the juvenile and adult dataframes together because there may be some juvenile taxa that weren't observed as adults or juveniles
#This code identifies which transects adult and juvenile colonies were recorded at and then converts NAs to 0s if needed
ssss<-subset(data.gen,GENUS_CODE=="SSSS")
ssss$Ad_pres<-ifelse(is.na(ssss$AdColCount),"0","-1")
ssss$Juv_pres<-ifelse(is.na(ssss$JuvColCount),"0","-1")
head(ssss)

ssss<-subset(ssss,select = c(SITE,SITEVISITID,TRANSECT,Ad_pres,Juv_pres,TRANSECTAREA_ad,TRANSECTAREA_j)) 
head(ssss)

data.gen<-left_join(subset(data.gen,select = -c(TRANSECTAREA_ad,TRANSECTAREA_j)),ssss) #use transect area from ssss because transectareas for some taxa were NA after merging adults and juvs

data.gen$JuvColCount[is.na(data.gen$JuvColCount)]<-0;data.gen$JuvColDen[is.na(data.gen$JuvColDen)]<-0
data.gen$AdColCount[is.na(data.gen$AdColCount) & data.gen$Ad_pres==-1]<-0;data.gen$AdColDen[is.na(data.gen$AdColDen) & data.gen$Ad_pres==-1]<-0
data.gen$Ad5_10ColCount[is.na(data.gen$Ad5_10ColCount)]<-0;data.gen$Ad5_10ColDen[is.na(data.gen$Ad5_10ColDen)]<-0


#Remove data from transects with less than 5m surveyed for adults and 1m for juvs.
data.gen$TRANSECTAREA_ad<-ifelse(data.gen$TRANSECTAREA_ad<5,NA,data.gen$TRANSECTAREA_ad);data.gen[data.gen$TRANSECTAREA_ad<5,]
data.gen$TRANSECTAREA_j<-ifelse(data.gen$TRANSECTAREA_j<1,NA,data.gen$TRANSECTAREA_j);data.gen[data.gen$TRANSECTAREA_j<1,]

#Calculate density of colonies <10cm
data.gen$Juv10ColDen<-data.gen$Ad5_10ColDen+data.gen$JuvColDen


#Site-level data
site.data.gen2<-subset(data.gen,TRANSECT==1) #drop transect 2 from 2015 data

head(site.data.gen2)


# GENERATE SUMMARY METRICS at the transect-leveL BY TAXONCODE--------------------------------------------------
#Calc_ColDen_Transect
acd.tax<-Calc_ColDen_Transect(data = awd,grouping_field = "TAXONCODE");colnames(acd.tax)[colnames(acd.tax)=="ColCount"]<-"AdColCount";colnames(acd.tax)[colnames(acd.tax)=="ColDen"]<-"AdColDen";colnames(acd.tax)[colnames(acd.tax)=="TRANSECTAREA"]<-"TRANSECTAREA_ad"# calculate density at genus level as well as total
jcd.tax<-Calc_ColDen_Transect(jwd,"TAXONCODE"); colnames(jcd.tax)[colnames(jcd.tax)=="ColCount"]<-"JuvColCount";colnames(jcd.tax)[colnames(jcd.tax)=="ColDen"]<-"JuvColDen";colnames(jcd.tax)[colnames(jcd.tax)=="TRANSECTAREA"]<-"TRANSECTAREA_j"

#Calc_ColMetric_Transect
cl.tax<-Calc_ColMetric_Transect(data = awd,grouping_field = "TAXONCODE",pool_fields = "COLONYLENGTH"); colnames(cl.tax)[colnames(cl.tax)=="Ave.y"]<-"Ave.cl" #Average % old dead
od.tax<-Calc_ColMetric_Transect(data = awd,grouping_field = "TAXONCODE",pool_fields = "OLDDEAD"); colnames(od.tax)[colnames(od.tax)=="Ave.y"]<-"Ave.od" #Average % old dead
rd.tax<-Calc_ColMetric_Transect(data = awd,grouping_field = "TAXONCODE",pool_fields = c("RDEXTENT1", "RDEXTENT2","RDEXTENT3")); colnames(rd.tax)[colnames(rd.tax)=="Ave.y"]<-"Ave.rd" #Average % recent dead

#ADD CODE TO CHANGE TRANSECT NUMBERS FOR JUVENILES
jcd.tax$TRANSECT[jcd.tax$TRANSECT==3]<-1
jcd.tax$TRANSECT[jcd.tax$TRANSECT==4]<-2

#Remove METHOD from dataframes before merging
acd.tax<-subset(acd.tax,select=-c(METHOD))
jcd.tax<-subset(jcd.tax,select=-c(METHOD))
cl.tax<-subset(cl.tax,select=-c(METHOD))
od.tax<-subset(od.tax,select=-c(METHOD))
rd.tax<-subset(rd.tax,select=-c(METHOD))


#Merge density and partial moratlity data together.You will need to replace the DUMMY field with the one you want
MyMerge <- function(x, y){
  df <- merge(x, y, by= c("SITE","SITEVISITID","TRANSECT","TAXONCODE"), all.x= TRUE, all.y= TRUE)
  return(df)
}
data.tax<-Reduce(MyMerge, list(acd.tax,jcd.tax,cl.tax,od.tax,rd.tax));


#Add METHOD back in
data.tax$METHOD<-"DIVER"

head(data.tax)

#Change NaN to NA
is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))

data.tax[is.nan(data.tax)] <- NA

#There will be some NAs when you merge the juvenile and adult dataframes together because there may be some juvenile taxa that weren't observed as adults or juveniles
#This code identifies which transects adult and juvenile colonies were recorded at and then converts NAs to 0s if needed
ssss<-subset(data.tax,TAXONCODE=="SSSS")
ssss$Ad_pres<-ifelse(is.na(ssss$AdColCount),"0","-1")
ssss$Juv_pres<-ifelse(is.na(ssss$JuvColCount),"0","-1")
head(ssss)

ssss<-subset(ssss,select = c(SITE,SITEVISITID,TRANSECT,Ad_pres,Juv_pres,TRANSECTAREA_ad,TRANSECTAREA_j)) 
head(ssss)

data.tax<-left_join(subset(data.tax,select = -c(TRANSECTAREA_ad,TRANSECTAREA_j)),ssss) #use transect area from ssss because transectareas for some taxa were NA after merging adults and juvs

data.tax$JuvColCount[is.na(data.tax$JuvColCount) & data.tax$Juv_pres==-1]<-0;data.tax$JuvColDen[is.na(data.tax$JuvColDen) & data.tax$Juv_pres==-1]<-0
data.tax$AdColCount[is.na(data.tax$AdColCount) & data.tax$Ad_pres==-1]<-0;data.tax$AdColDen[is.na(data.tax$AdColDen) & data.tax$Ad_pres==-1]<-0


#Remove data from transects with less than 5m surveyed for adults and 1m for juvs.
data.tax$TRANSECTAREA_ad<-ifelse(data.tax$TRANSECTAREA_ad<5,NA,data.tax$TRANSECTAREA_ad);data.tax[data.tax$TRANSECTAREA_ad<5,]
data.tax$TRANSECTAREA_j<-ifelse(data.tax$TRANSECTAREA_j<1,NA,data.tax$TRANSECTAREA_j);data.tax[data.tax$TRANSECTAREA_j<1,]

#Site-level data
site.data.tax2<-subset(data.tax,TRANSECT==1) #drop transect 2 from 2015 data

head(site.data.tax2)

#Merge Site Master with demographic data
sm<-read.csv("T:/Benthic/Projects/Swains 2023 Benthic Analysis/SURVEY_MASTER_w2023benthic.csv", stringsAsFactors=FALSE)
site.data.gen2<-left_join(site.data.gen2,sm[,c("OBS_YEAR","ISLAND","SEC_NAME","SITEVISITID","SITE","DEPTH_BIN","new_MIN_DEPTH_M","new_MAX_DEPTH_M")]) 
site.data.tax2<-left_join(site.data.tax2,sm[,c("OBS_YEAR","ISLAND","SEC_NAME","SITEVISITID","SITE","DEPTH_BIN","new_MIN_DEPTH_M","new_MAX_DEPTH_M")]) 

#Write files
write.csv(file="T:/Benthic/Projects/Swains 2023 Benthic Analysis/Data/Swains_sitedata_GENUS.csv",site.data.gen2)
write.csv(file="T:/Benthic/Projects/Swains 2023 Benthic Analysis/Data/Swains_sitedata_TAXONCODE.csv",site.data.tax2)

