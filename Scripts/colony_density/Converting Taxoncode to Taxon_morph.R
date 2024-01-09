library(dplyr)

data <- read.csv("CoralBelt_Adults_raw_CLEANED.csv")

data$TAXONCODE_2 <- NA

##### POSP ######
# merge Foliose, Laminar, Mounding, Plating, Knobby, Columnar 
data[data$GENUS_CODE == "POSP" & 
       data$MORPHOLOGY %in% c("Columnar","Foliose","Knobby",
                              "Laminar","MD","Mounding","PL","Plating"),]$TAXONCODE_2 <- "POSP_MD"
# merge Encrusting flat + mounding 
data[data$GENUS_CODE == "POSP" & 
       data$MORPHOLOGY %in% c("EM","Encrusting (flat)","Encrusting (mounding)"),]$TAXONCODE_2 <- "POSP_EM"
# drop 4 branching individuals 
data <- data[- which(data$GENUS_CODE == "POSP" & 
               data$MORPHOLOGY %in% c("BR","Branching")),]


##### MOSP ######
# merge mounding, EF, & EM
data[data$GENUS_CODE == "MOSP" & 
       data$MORPHOLOGY %in% c("EM","Encrusting (flat)","Encrusting (mounding)","MD","Mounding"),]$TAXONCODE_2 <- "MOSP_EM"
# merge TB, FO, LM & PL
data[data$GENUS_CODE == "MOSP" & 
       data$MORPHOLOGY %in% c("FO","Foliose","Laminar","PL","Plating","Tabulate"),]$TAXONCODE_2 <- "MOSP_FO"
# drop 7 branching, 1 knobby, and 2 unidentified morphs
data <- data[-which(data$GENUS_CODE == "MOSP" & 
               data$MORPHOLOGY %in% c("BR","Branching","KN","_ ")),]


##### POCS ######
# all POCS are branching
# PMVC (PDAM, PVER, PMEA, PMVC, PDAN) 
data[data$SPCODE %in% c("PMVC","PDAM","PVER","PMEA","PDAN","POCS"),]$TAXONCODE_2 <- "PMVC"
# PGWC (PEYD, PWOO, PGWC)
data[data$SPCODE %in% c("PGWC","PEYD","PWOO"),]$TAXONCODE_2 <- "PGWC"


##### PSSP ######
# roll all up to genus-morphology given confusion of species during 23
# merge knobby, branching, and columnar 
data[data$GENUS_CODE == "PSSP" & 
       data$MORPHOLOGY %in% c("BR","CO","KN","Knobby"),]$TAXONCODE_2 <- "PSSP_KN"
# merge encrusting all together, with mounding 
data[data$GENUS_CODE == "PSSP" & 
       data$MORPHOLOGY %in% c("EM","Encrusting (columnar)","Encrusting (mounding)","Mounding"),]$TAXONCODE_2 <- "PSSP_EM"
# drop 1 foliose 
data <- data[- which(data$GENUS_CODE == "PSSP" & 
               data$MORPHOLOGY %in% c("FO")),]


##### PAVS ######
# keep PCHI, PDUE, PMAL, PVAR as species 
data[data$SPCODE == "PCHI",]$TAXONCODE_2 <- "PCHI"
data[data$SPCODE %in% c("PDUE","PDIF"),]$TAXONCODE_2 <- "PDUE"
data[data$SPCODE == "PMAL",]$TAXONCODE_2 <- "PMAL"
data[data$SPCODE == "PVAR",]$TAXONCODE_2 <- "PVAR"
# change PAVS to encrusting morphology
data[data$SPCODE == "PAVS",]$TAXONCODE_2 <- "PAVS_EM"


###### REMOVE ALL FUSP ######
data <- data[- which(data$GENUS_CODE == "FUSP"),]

###### LESP ######
# leave to species level
data[data$SPCODE == "LINC",]$TAXONCODE_2 <- "LINC"
data[data$SPCODE == "LMYC",]$TAXONCODE_2 <- "LMYC"
# drop 1 LESP
data <- data[- which(data$SPCODE == "LESP"),]


###### COSP ######
# all at genus level, encrusting 
data[data$GENUS_CODE == "COSP",]$TAXONCODE_2 <- "COSP_EM"


###### STYS ######
# keep as is (all branching)
data[data$GENUS_CODE == "STYS",]$TAXONCODE_2 <- "SPIS"


###### HESP ######
# keep as is 
data[data$GENUS_CODE == "HESP",]$TAXONCODE_2 <- "HCOE"


###### FASP ######
# keep as is (genus code)
data[data$GENUS_CODE == "FASP",]$TAXONCODE_2 <- "FSTE"


###### ACSP ######
# Roll all up to genus code (ACSP branching)
data[data$GENUS_CODE == "ACSP",]$TAXONCODE_2 <- "ACSP_BR"


###### CYSP ######
# Not using for this analysis, but all need to be FUSP
data <- data[- which(data$GENUS_CODE == "CYSP"),]


###### PLSP ######
# keep as is (genus code)
data[data$GENUS_CODE == "PLSP",]$TAXONCODE_2 <- "PLSP"


###### ASTS ######
# keep as is (genus code)
data[data$GENUS_CODE == "ASTS",]$TAXONCODE_2 <- "ASTS"


###### THE REST ######
data$TAXONCODE_2 <- ifelse(is.na(data$TAXONCODE_2), data$SPCODE, data$TAXONCODE_2)
# removing UNKN
data <- data[- which(data$TAXONCODE_2 == "UNKN"),]

# Double check 
unique(data$TAXONCODE_2)


write.csv(data, "CoralBelt_Adults_raw_CLEANED_Swains_GenusMorph.csv")



