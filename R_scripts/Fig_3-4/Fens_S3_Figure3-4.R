#' ---
#' eDNA matching sites Sp Richness"
#' Nathan P Griffiths
#' July 2020
#' -------------------------------------
#-------------------------------SETUP & Raw Data Input-------------------------------------------#

## Clear Memory
rm(list=ls())

## Check Working Directory 
getwd()

## To ensure reproducibility, print R.version details.
sessionInfo()
#R version 3.6.3 (2020-02-29)
#Platform: x86_64-w64-mingw32/x64 (64-bit)
#Running under: Windows 10 x64 (build 18362)

## Load catch data 
catch <- read.csv("eel2017_catch_DNAlayout.csv", header=TRUE)
head(catch)

## Load eDNA data 
edna <- read.csv("eel2017_sorted.csv", header=TRUE)
head(edna)


#### Organise Data ####

library(vegan)
library(dplyr)

#set edna dimentions as r
r=dim(edna)

#Remove assigned column
edna=edna[1:(r[2]-1)]

#Make new dataframe with first column as row names
edna1 <- data.frame(edna[,-1], row.names=edna[,1])
catch1 <- data.frame(catch[,-1], row.names=catch[,1])

#Make edna proportional
edna1$assigned=rowSums(edna1)
edna1=edna1/edna1$assigned
r=dim(edna1)
edna1=edna1[1:(r[2]-1)]

#Make catch proportional
catch1$assigned=rowSums(catch1)
catch1=catch1/catch1$assigned
r=dim(catch1)
catch1=catch1[1:(r[2]-1)]

##Compare Species Richness Visually 
edna1$richness = rowSums(edna1>0)
plot(edna1$richness,ylim=c(0,20))

catch1$richness = rowSums(catch1>0)
plot(catch1$richness,ylim=c(0,20))

#### Make Dataframe For Plot 
y = data.frame(edna$X, catch1$richness, edna1$richness)
colnames(y) <- c("Site","Catch","eDNA")

library(reshape2)

det = melt(y, id=(1))

## ggplot stuff ## 
det$Site<-as.factor(det$Site)
library(ggplot2)

##Compare Species Richness Visually 
p <- ggplot(data=det, aes(x=Site, y=value, group=variable)) +
  geom_line(size=1, aes(color=variable)) +
  geom_point(size=2, aes(color=variable)) + 
  scale_color_manual(values=c("darkred", "darkblue"))

p + theme_classic() + 
  theme(legend.position = "bottom") +
  xlab("Site Number") + 
  ylab ("Species Richness") +
  theme(legend.title = element_blank())

## Scatter plot ## 

## save det as df for restructure to plot as scatter plot ## 

## Save dataframe as .csv 
write.csv(det, file = 'eel2017_scatter.csv')

## Modify "eel2017_scatter" so eDNA & Catch are as seperate columns. 

## Read in modified data ## 
sp <- read.csv("SpR.csv", header=TRUE)
head(sp)

library(ggrepel)
## Make Scatterplot 
z <- ggplot (sp, aes(x=Catch, y=eDNA)) + 
  geom_point(size=2, shape=19, color="Black") +
  xlab("Species Richness (Catch)") +
  ylab("Species Richness (eDNA)") +
  theme_bw() +
  geom_segment(aes(x = 0, y = 0, xend = 20, yend = 20), size=0.5, colour="darkred", linetype="longdash") +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0))

## Save Scatterplot
png('scattersite.png',height=1500,width=1500,res=350)
par(mar=c(1,1,1,1),xpd=T)
z
dev.off()

## Text Label points
zln = ggplot (sp, aes(x=Catch, y=eDNA)) + 
  xlab("Species Richness (Catch)") +
  ylab("Species Richness (eDNA)") +
  theme_bw() +
  geom_segment(aes(x = 0, y = 0, xend = 20, yend = 20), size=0.5, colour="darkred", linetype="longdash") +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) + 
  geom_label_repel(aes(label = Site),
                   size = 3,
                   box.padding   = 0, 
                   point.padding = 0)

png('scattersiteL.png',height=1500,width=1500,res=350)
par(mar=c(1,1,1,1),xpd=T)
zln
dev.off()

## Make Species Level Scatterplot ## 

spscatch <- catch
spsedna <- edna

#Make new dataframe with first column as row names
ea <- data.frame(spscatch[,-1], row.names=spscatch[,1])
ed <- data.frame(spsedna[,-1], row.names=spsedna[,1])

#Transpose Data with site as column names
ea1 <- data.frame(t(ea)) 
ed1 <- data.frame(t(ed))

## Make column for species count 
ea1$sc = rowSums(ea1>0)
ed1$sc = rowSums(ed1>0)

## Save dataframes as .csv 
write.csv(ea1, file = 'ea1.csv')
write.csv(ed1, file = 'ed1.csv')

## Modify "ea1, ed1" so eDNA & Catch are merged as seperate columns. 

## Read in modified data ## 
so <- read.csv("SpO.csv", header=TRUE)
head(so)

## Make Scatterplot 
o <- ggplot (so, aes(x=Catch, y=eDNA)) + 
  geom_point(size=2, shape=19, color="Black") +
  xlab("Species Site Occupancy (Catch)") +
  ylab("Species Site Occupancy (eDNA)") +
  theme_bw() +
  geom_segment(aes(x = 0, y = 0, xend = 20, yend = 20), size=0.5, colour="darkred", linetype="longdash") +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) 
  #geom_jitter(height=0, width=0.5 , size=2, shape=19, color="Black")

png('scatteroccupancy.png',height=1500,width=1500,res=350)
par(mar=c(1,1,1,1),xpd=T)
o
dev.off()

## Text Labels >10
m = ggplot (so, aes(x=Catch, y=eDNA)) + 
  geom_point(size=2, shape=19, color="Black") +
  xlab("Species Site Occupancy (Catch)") +
  ylab("Species Site Occupancy (eDNA)") +
  theme_bw() +
  geom_segment(aes(x = 0, y = 0, xend = 20, yend = 20), size=0.5, colour="darkred", linetype="longdash") +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  geom_label_repel(aes(label=ifelse(eDNA>10,as.character(abbreviate(Species, minlength = 4)),'')),
                   size = 3,
                   box.padding   = 0, 
                   point.padding = 0)

## Label only eel 
es = ggplot (so, aes(x=Catch, y=eDNA)) + 
  geom_point(size=2, shape=19, color="Black") +
  xlab("Species Site Occupancy (Catch)") +
  ylab("Species Site Occupancy (eDNA)") +
  theme_bw() +
  geom_segment(aes(x = 0, y = 0, xend = 20, yend = 20), size=0.5, colour="darkred", linetype="longdash") +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  geom_label_repel(aes(label=ifelse(eDNA>13 & Catch<5,as.character(abbreviate(Species, minlength = 4)),'')),
                   size = 3,
                   box.padding   = 0, 
                   point.padding = 0)

png('scattereel.png',height=1500,width=1500,res=350)
par(mar=c(1,1,1,1),xpd=T)
es
dev.off()


####################### Statistical Tests #######################

## Richness ## #sp 

##Check for correlation##
#Test for Normality
qqnorm(sp$Catch)
qqnorm(sp$eDNA)
hist(sp$Catch)
hist(sp$eDNA)
shapiro.test(sp$Catch)
shapiro.test(sp$eDNA)
#Data are normal - Pearson 
cor.test(~ Catch + eDNA, data = sp,
         method = "pearson",
         conf.level = 0.95)
#Significant Correlation
#Plot
library("ggpubr")
cr <- ggscatter(sp, x= "Catch", y ="eDNA",
          add = "reg.line", conf.int = FALSE,
          cor.coef = TRUE, cor.method = "pearson",  cor.coef.coord = c(1, 17),
          xlab = "Species Richness (Catch)",
          ylab = "Species Richness (eDNA)")
png('Correlation_Richness.png',height=1500,width=1500,res=350)
par(mar=c(1,1,1,1),xpd=T)
cr
dev.off()

#If data were not normal - Spearman 
cor.test( ~ Catch + eDNA,
          data=sp,
          method = "spearman",
          continuity = FALSE,
          conf.level = 0.95,
          exact=FALSE)
#Significant Correlation

########################

#calculate the mean species richness for Catch and eDNA. 
mean(sp$Catch)
mean(sp$eDNA)

# calculate standard error of the mean
sd(sp$Catch)/sqrt(length(sp$Catch))
sd(sp$eDNA)/sqrt(length(sp$eDNA))

#Melt Test 
mtest <- sp[,-1] 
mtest <- melt(mtest)

#table of means 
table(mtest$variable)
tapply(mtest$value, mtest$variable, mean)

## Paired t-test 
G1 <- subset(mtest, variable=='Catch')
G2 <- subset(mtest, variable=='eDNA')
t.test(G1$value, G2$value, paired=T)

#Differences between paired values are normal.
qqnorm(sp$Catch - sp$eDNA)
shapiro.test(sp$Catch - sp$eDNA)
hist(sp$Catch - sp$eDNA)

#If data were non normal#
#Paired test 
wilcox.test(sp$Catch, sp$eDNA, paired=TRUE)

#ggplot ## Boxplot ## 
ta <- melt(sp, id = c("Site"))

library(ggsignif)

br <- ggplot(ta, aes(x=variable, y=value, fill=variable))+
  stat_boxplot(geom = "errorbar", width = 0.2)+
  geom_boxplot()+
  geom_jitter(width=0.15, height=0.1, size=1)+
  ylab("Species Richness")+
  xlab("Survey Method")+
  scale_x_discrete() +
  theme_classic() +
  geom_signif(y_position = 20,
              xmin = 1, xmax = 2, 
              annotation = "p < 0.001", tip_length = 0.02) +   
scale_fill_brewer(palette="Greys") +
  theme(legend.position = "none")

png('boxplotrichness.png',height=1500,width=1500,res=350)
par(mar=c(1,1,1,1),xpd=T)
br
dev.off()

## Merge Richness Plot 
library(gridExtra)
f3 <- grid.arrange(z, br, cr, ncol = 3, nrow = 1)

library(cowplot)
fin3 <- ggdraw(f3) + draw_plot_label(c("(a)","(b)", "(c)"), 
                                   c(0.031, 0.37, 0.71), 
                                   c(0.99, 0.99, 0.99), 
                                   fontface = "plain",size = 11)

##################################################

## Site Occupancy ## so

##Check for correlation##
#Test for Normality
qqnorm(so$Catch)
qqnorm(so$eDNA)
hist(so$Catch)
hist(so$eDNA)
shapiro.test(so$Catch)
shapiro.test(so$eDNA)
#Data are not normal - Spearman 
cor.test( ~ Catch + eDNA,
          data=so,
          method = "spearman",
          continuity = FALSE,
          conf.level = 0.95,
          exact=FALSE)
#Significant Correlation
#Plot
co <- ggscatter(so, x= "Catch", y ="eDNA",
                add = "loess", conf.int = FALSE,
                cor.coef = TRUE, cor.method = "spearman", cor.coef.coord = c(1, 17),
                xlab = "Species Site Occupancy (Catch)",
                ylab = "Species Site Occupancy (eDNA)")
png('Correlation_Occupancy.png',height=1500,width=1500,res=350)
par(mar=c(1,1,1,1),xpd=T)
co
dev.off()

##If treated as normal - Pearson
cor.test(~ Catch + eDNA, data = so,
         method = "pearson",
         conf.level = 0.95)
#Significant Correlation

############################################

#calculate the mean Site Occupancy for Catch and eDNA. 
mean(so$Catch)
mean(so$eDNA)

# calculate standard error of the mean
sd(so$Catch)/sqrt(length(sp$Catch))
sd(so$eDNA)/sqrt(length(sp$eDNA))

#Melt Test 
mtest <- so[,-1] 
mtest <- mtest[c(2,1)]
mtest <- melt(mtest)

#table of means 
table(mtest$variable)
tapply(mtest$value, mtest$variable, mean)

## Paired t-test 
G1 <- subset(mtest, variable=='Catch')
G2 <- subset(mtest, variable=='eDNA')
t.test(G1$value, G2$value, paired=T)

#Differences between paired values look normal.
qqnorm(so$Catch - so$eDNA)
shapiro.test(so$Catch - so$eDNA)
hist(so$Catch - sp$eDNA)

#Data look to conform to normality, however shapiro.test suggests they may differ
#from a normal distribution. Paired wilcox.test gives same output. 
#Paired test 
wilcox.test(so$Catch, so$eDNA, paired=TRUE, exact = FALSE)

#ggplot ## Boxplot ## 
so <- so[c(1,3,2)]
ta <- melt(so, id = c("Species"))

library(ggsignif)


bo <- ggplot(ta, aes(x=variable, y=value, fill=variable))+
  stat_boxplot(geom = "errorbar", width = 0.2)+
  geom_boxplot()+
  geom_jitter(width=0.15, height=0.1, size=1)+
  ylab("Site Occupancy")+
  xlab("Survey Method")+
  theme_classic() +
  theme(legend.position = "none") +
  scale_x_discrete() +
  scale_fill_brewer(palette="Greys") +
geom_signif(y_position = 20,
            xmin = 1, xmax = 2, 
            annotation = "p < 0.01", tip_length = 0.02)  

png('boxplotoccupancy.png',height=1500,width=1500,res=350)
par(mar=c(1,1,1,1),xpd=T)
bo
dev.off()
#####################################################################

## Merge Occupancy 

f4 <- grid.arrange(es, bo, co, ncol = 3, nrow = 1)

fin4 <- ggdraw(f4) + draw_plot_label(c("(a)","(b)", "(c)"), 
                                     c(0.031, 0.37, 0.71), 
                                     c(0.99, 0.99, 0.99), 
                                     fontface = "plain",size = 11)

######################################################################

## Figure 3
c3 = cr +
  geom_segment(aes(x = 0, y = 0, xend = 20, yend = 20), size=0.5, colour="darkred", linetype="longdash") + theme_bw() +
  scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0))

library(gridExtra)
f3 <- grid.arrange(br, c3, ncol = 2, nrow = 1)

library(cowplot)
fin3 <- ggdraw(f3) + draw_plot_label(c("(a)","(b)"), 
                                     c(0.05, 0.55), 
                                     c(0.98, 0.98), 
                                     fontface = "plain",size = 11)

png('merge_fig3.png',height=1500,width=3000,res=350)
par(mar=c(1,1,1,1),xpd=T)
fin3
dev.off()

## Figure 4
c4 = co +
  geom_segment(aes(x = 0, y = 0, xend = 20, yend = 20), size=0.5, colour="darkred", linetype="longdash") + theme_bw() + 
  geom_label_repel(aes(label=ifelse(eDNA>13 & Catch<5,as.character(abbreviate(Species, minlength = 4)),'')),
                   size = 3,
                   box.padding   = 0, 
                   point.padding = 0)

f4 <- grid.arrange(bo, c4, ncol = 2, nrow = 1)

fin4 <- ggdraw(f4) + draw_plot_label(c("(a)","(b)"), 
                                     c(0.05, 0.55), 
                                     c(0.98, 0.98), 
                                     fontface = "plain",size = 11)

png('merge_fig4.png',height=1500,width=3000,res=350)
par(mar=c(1,1,1,1),xpd=T)
fin4
dev.off()

