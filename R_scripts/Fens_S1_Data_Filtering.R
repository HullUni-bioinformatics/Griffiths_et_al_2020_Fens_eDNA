## Clean up for raw metaBEAT output ##
#' Nathan P Griffiths - Eels in the fens 2017 
#' July 2020
#' -----------------------------------

##################### Prepare working environment #####################

## Clear memory
rm(list=ls())

## Check working directory
getwd()

## Set working directory if required 
#setwd("file path")

## To ensure reproducibility, print details about the version of R being
## used for analysis.
sessionInfo()

#R version 3.6.3 (2020-02-29)
#Platform: x86_64-w64-mingw32/x64 (64-bit)
#Running under: Windows 10 x64 (build 18362)

## Load raw data 
dat_raw <- read.csv("Fens_2017_Fish_0.98_Raw.csv", header=TRUE)
head(dat_raw)
#-------------------------------Raw Data Processing-------------------------------------------#

## Reset row names
rownames(dat_raw) <- NULL

## Remove 'X' and '.nc.blast' from Site Names 
colnames(dat_raw) <- gsub('X', '', colnames(dat_raw))
colnames(dat_raw) <- gsub('.nc.blast', '', colnames(dat_raw))


## Rename first & last column
r=dim(dat_raw)
colnames(dat_raw)[1] <- "Assignment"
colnames(dat_raw)[r[2]] <- "Tax"

## Transpose data & remove taxonomy column 
dat_t=data.frame(t(dat_raw[,!colnames(dat_raw) %in% 'Tax',][,-1]))
colnames(dat_t)=dat_raw[,1]
dat_t
str(dat_t)

## Add assigned reads column 
dat_t -> me
r=dim(me)
r[2]
me$assigned=rowSums(me[1:(r[2]-1)])
head(me)

## Add total reads column
me$total=me$assigned+me$unassigned
head(me)

## Remove unassigned reads column
me=me[,!colnames(me)%in%'unassigned',]
head(me)
dim(me)

## Remove columns with 0 reads 
me=me[, colSums(me != 0) > 0]
head(me)
dim(me)

## Save dataframe as .csv
write.csv(me, file = "Fens_uncleaned.csv",row.names=TRUE)

##### Apply a threshold to remove reads <0.1% to account for contaminants/false reads #####

## Create 'dat_p': a proportion reads dataframe
dat_p=me/me$total
dat_p
dim(dat_p)

## Make test dataframe 
test1 <- me

## Apply 0.1% threshold from dat to test1
test1[dat_p<0.001]=0

##### Finalise dataframe #####

#' 0.1% threshold applied 

## Rename test1 
dat <- test1

## Remove control Sp. 
dat=dat[,!colnames(dat)%in%'Maylandia_zebra',]
dim(dat)

## Remove columns with no reads after threshold 
dat=dat[, colSums(dat != 0) > 0]
dim(dat)

## Recalculate assigned reads - They have changed due to appication of threshold
dat=dat[,!colnames(dat)%in%'assigned',]
r=dim(dat)
dat$assigned=rowSums(dat[1:(r[2]-1)])
head(dat)
dim(dat)

## Save dataframe 
write.csv(dat, file = "Fens_threshold.csv",row.names=TRUE)

#------------------------------- Refine Dataset -------------------------------------------#

## Remove family level assignments (non-specific families)
dat=dat[,!colnames(dat)%in%'Coregonus',]
dat=dat[,!colnames(dat)%in%'Cyprinidae',]
dat=dat[,!colnames(dat)%in%'Gasterosteidae',]
dat=dat[,!colnames(dat)%in%'Salmonidae',]

#' Now, correct species names:
#' 
#' Lampetra_fluviatilis = Lampetra
#' Cobitidae = Cobitis_taenia
#' Percidae = Perca_fluviatilis
library(plyr)
test <- rename(dat, c("Lampetra_fluviatilis"="Lampetra", "Percidae"="Perca_fluviatilis"))

## Now, merge any tax with same names before moving forward
str(test)
test$testt <- test$Cobitis_taenia + test$Cobitidae
str(test)

## Delete columns 
test=test[,!colnames(test)%in%'Cobitis_taenia',]
test=test[,!colnames(test)%in%'Cobitidae',]

## Rename sum of two columns 
library(dplyr)
str(test)
t <- rename(test, Cobitis_taenia = testt)
str(t)

dat <- t 

#Create row of occupancy for each species 
dat["Occupancy" ,] <- colSums(dat>0)

#Remove columns with occupancy =x, to account for environmental contaminants and anomalies
dim(dat) #Get dimentions of df
t=dat[107,] >1 #Make t = data to keep 
t
p <- dat #Make p = df

#t = p[,!(p[107,]==1)] #Remove all species with only 1 occurance in df
#t = p[,!(p[107,]<2)] #Remove all species with <x occurances in df 
#' We keep all species for further analysis, but flag Carassius_auratus, Leuciscus_idus, & Oncorhynchus_mykiss
#' which only occur in one sample as possible environmental contamination. 
t = p 

## Remove occupancy row 
datt <- t[-107, ]

## Remove total & assigned columns
datt=datt[,!colnames(datt)%in%'total',]
datt=datt[,!colnames(datt)%in%'assigned',]

## Remove columns that = 0
datt=datt[, colSums(datt != 0) > 0]
dim(datt)

## Save dataframe 
write.csv(datt, file = "Fens_2017_cleaned.csv",row.names=TRUE)

