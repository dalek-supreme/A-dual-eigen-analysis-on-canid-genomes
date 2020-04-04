library(ape)
snp.tree <- read.tree("snphylo.output.ml.tree")
#plot(snp.tree)
library(tidytree)

# to map
library(phytools)
demographic.data <- read.csv('~/demographic.csv')
demographic.data.valid <- demographic.data[!(is.na(demographic.data$Latitude)|is.na(demographic.data$Longitude)),]

lat.long <- data.frame(
    lat=demographic.data.valid$Latitude,
    long=demographic.data.valid$Longitude
)

row.names(lat.long) <- demographic.data.valid$SampleID

invalid.id <- as.character(demographic.data[is.na(demographic.data$Latitude)|is.na(demographic.data$Longitude),]$SampleID)
#valid.id <- as.character(demographic.data[!(is.na(demographic.data$Latitude)|is.na(demographic.data$Longitude)),]$SampleID)

obj<-phylo.to.map(drop.tip(snp.tree,invalid.id),lat.long,plot=FALSE)

plot(obj,type="phylogram",asp=1)