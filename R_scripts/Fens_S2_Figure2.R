#' ---
#' eDNA matching sites Sp composition test
#' Nathan P Griffiths
#' July 2020
#' -------------------------------------

#'The output from fens_S1 has been subset to include only the 17 sites which correspond 
#'to traditional surveys carried out that same year. The average eDNA reads for each 
#'site have been calculated from the 5 samples, and are compared with catch 
#'methods here. This was done outside of R to enable manual checks of the data during this
#'process. The eDNA_Matching_Sites.csv file is input here. 

#-------------------------------SETUP & Raw Data Input-------------------------------------------#

## Clear Memory
rm(list=ls())

## Check Working Directory 
getwd()

## Load data 
raw <- read.csv("eDNA_Matching_Sites_AllSp.csv", header=TRUE)
head(raw)

## To ensure reproducibility, print R.version details.
sessionInfo()
#R version 3.6.3 (2020-02-29)
#Platform: x86_64-w64-mingw32/x64 (64-bit)
#Running under: Windows 10 x64 (build 18362)

#-------------------------------Raw Data Processing-------------------------------------------#

## Reset row names
rownames(raw) <- NULL

## Rename first column
colnames(raw)[1] <- "Site"

## Remove columns with 0 reads 
raw=raw[, colSums(raw != 0) > 0]
head(raw)
dim(raw)

## Save dataframe as .csv
write.csv(raw, file = "eel2017_Sp.csv",row.names=FALSE)

###FRESH START TO ORGANISE FOR PLOTS###
## Read.csv
plt=read.csv('eel2017_Sp.csv',header=TRUE, row.names=1)
head(plt)

## Now remove the underscore from species names
colnames(plt) <- gsub('_', ' ', colnames(plt))

## Sort column names alphabetically
plt <- plt[order(names(plt))]
names(plt)

## Add back assigned column
r=dim(plt)
plt$assigned=rowSums(plt[1:(r[2])])
head(plt)

##ggplot stuff## 
## Load Library Reshape to 'melt' data and enable data to fit ggplot parameters 
library(reshape2)
library(tidyverse)

## Save dataframe as .csv 
write.csv(plt, file = 'eel2017_sorted.csv')



############################################################################
## Bubble Plots ## 

## Make Plot - all eDNA Matching sites ##
# Libraries
library(dplyr)
library(hrbrthemes)
library(viridis)
library(ggplot2)


## Make proportional reads eDNA dataset 
dim(plt)
site = rownames(plt)
det <- plt
r=dim(det)
det=det[1:(r[2]-1)]
det$assigned=rowSums(det)
det=det/det$assigned
r=dim(det)
det=det[1:(r[2]-1)]
det = cbind(site, det) 
dim(det)
det = melt(det, id=(1))
View(det)

## Save dataframe as .csv 
write.csv(det, file = 'eel2017_melted.csv')


## Make eDNA Bubble Plot 

png('eDNA_Bubb.png',height=1800,width=2700,res=300)
par(mar=c(1,1,1,1),xpd=T)

Set1 <- subset(det,det$value>0)
str(Set1)
Set1$site <- factor(Set1$site, levels = site)
Set1$site
bubb1 <- ggplot(Set1, aes(x=site, y=variable, size=value*100)) +
  geom_point(alpha= .8, shape=21, fill="darkgreen", color = "black") +
  scale_size(range = c(1, 12), name="% Reads") +
  theme_bw() +
  theme(legend.position = "right") +
  xlab("Site Number") +
  ylab("") +
  #xlab("Site/Sample") +
  #ylab("Species") +
  theme(axis.text.y = element_text(face = "italic", size = 12)) +
  theme(axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 14)) +
  scale_x_discrete(position = "top") +
  scale_y_discrete(limits = rev(levels(Set1$variable)))
bubb1

dev.off()


## Make Plot - all Environment Agency Matching sites ##
## Load data 
raw <- read.csv("EA_DATA_MATCH.csv", header=TRUE)
head(raw)

dim(raw)
det <- raw
r=dim(det)
View(det)

Pet1 <- subset(det,det$ALL_RUNS>0)
str(Pet1)
Pet1$SITE_N <- as.factor(Pet1$SITE_N)


######################## EA % CATCH #######################

## Unmelt Data 
TEST <- dcast(Pet1, SITE_N ~ LATIN_NAME, value.var = "ALL_RUNS")

## Remove na's
TEST[is.na(TEST)] <- 0 

## Reset row names
rownames(TEST) <- NULL

## Remove columns with 0 reads 
TEST=TEST[, colSums(TEST != 0) > 0]
head(TEST)
dim(TEST)

rownames(TEST) <- TEST$SITE_N

r=dim(TEST)
r

#Remove site column
t=TEST[,-1]

## Save dataframe as .csv 
write.csv(t, file = 'eel2017_catch_DNAlayout.csv')

## Make proportional catch dataset 
dim(t)
site = rownames(t)
det <- t
r=dim(det)
det$assigned=rowSums(det)
det=det/det$assigned
r=dim(det)
det=det[1:(r[2]-1)]
det = cbind(site, det) 
dim(det)
det = melt(det, id=(1))
View(det)

## Remove 0's 
Pet1 <- subset(det,det$value>0)
str(Pet1)

## Save dataframe as .csv 
write.csv(Pet1, file = 'eel2017_catch_melted.csv')


##### Manually add back in survey methods and order site column #####


## Load sorted data 
Pet1 <- read.csv("2017eelcatch_sorted.csv", header=TRUE)
str(Pet1)
Pet1$site <- as.factor(Pet1$site)


## Make Catch Bubble Plot

png('EA_bubb.png',height=1800,width=2900,res=300)
par(mar=c(1,1,1,1),xpd=T)


bubb2 <- ggplot(Pet1, aes(x=site, y=variable, size=value*100)) +
  geom_point(alpha= .8, shape=21, fill="darkblue", color="black") +
  scale_size(range = c(1, 12), name="% Catch") +
  theme_bw() +
  theme(legend.position = "right") +
  xlab("Site Number") +
  ylab("") +
  #xlab("Site/Sample") +
  #ylab("Species") +
  theme(axis.text.y = element_text(face = "italic", size = 12)) +
  theme(axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 14)) +
  scale_x_discrete(position = "top")  +
  scale_y_discrete(limits = rev(levels(Pet1$variable)))
bubb2

dev.off()


## Merge Plots 
library(gridExtra)
grid.arrange(bubb1, bubb2, ncol=2)

## Standardise and index df's
e1=Set1
c1=Pet1[,-3]
e1$id <- rep(1, nrow(e1))
c1$id <- rep(2, nrow(c1))

## Merge
m1 <- rbind(e1, c1)
m1$id <- as.factor(m1$id)
m1$variable <- as.character(m1$variable)
str(m1$id)


bubb3 <- ggplot(m1, aes(x=site, y=fct_rev(variable), size=value*100, fill=id)) +
  geom_point(alpha= .8, shape=21, color="black") +
  scale_size(range = c(1, 14), name="% Reads/Catch") +
  scale_fill_manual(values = c("darkgreen", "darkblue")) +
  theme_bw() +
  theme(legend.position = "right") +
  xlab("Site Number") +
  ylab("") +
  theme(axis.text.y = element_text(face = "italic", size = 12)) +
  theme(axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 14)) +
  scale_x_discrete(position = "top")  +
  scale_y_discrete(limits = rev(levels(m1$variable)))
bubb3 

p <- bubb3 + guides(fill = FALSE)

fig2 = p + facet_grid(cols = vars(id),
               labeller = as_labeller(c("1"="eDNA", "2"="Catch"))) +
  theme(
    strip.text.x = element_text(
      size = 12, face = "bold"))


##SAVE FIGURE 2 
png('FIG2.png',height=2300,width=4200,res=300)
par(mar=c(1,1,1,1),xpd=T)
fig2
dev.off()



