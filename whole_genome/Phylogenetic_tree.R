library(ape)
snp.tree <- read.tree("snphylo.output.ml.tree")
#plot(snp.tree)
#library(tidytree)

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

obj<-phylo.to.map(drop.tip(snp.tree,invalid.id),lat.long,plot=FALSE)

plot(obj,type="phylogram",asp=1)


# python TreeCluster.py -i snphylo.output.ml.tree2 -t 0.4 > cluster_results 
cluster.result <- read.table('cluster_results',header=T)

cluster.result <- cluster.result[cluster.result$SequenceName!='Lcu2_Pasto',]
# remove invalid item


cluster.map <- data.frame(
    SampleID=cluster.result$SequenceName,
    Latitude=array(0,length(cluster.result$SequenceName)),
    Longitude=array(0,length(cluster.result$SequenceName)),
    ClusterNumber=cluster.result$ClusterNumber
)


for (i in rownames(cluster.map)) {
    demog.row <- demographic.data.valid[as.character(demographic.data.valid$SampleID)==as.character(cluster.map[i,]$SampleID),]
    cluster.map[i,]$Latitude <- demog.row$Latitude
    cluster.map[i,]$Longitude <- demog.row$Longitude
}

#plot on map
library(ggplot2)
library(ggmap)
library(sp)
library(maptools)
library(maps)


mapworld<-borders("world",colour = "gray50",fill="white")
mp<-ggplot()+mapworld+ylim(-60,90)
mp<-mp+geom_point(aes(x=as.numeric(cluster.map$Longitude),
                        y=as.numeric(cluster.map$Latitude),color=as.factor(cluster.map$ClusterNumber)))
# mid<-mean(unique(cluster.map$ClusterNumber))
# mp<-mp+scale_color_gradient2(midpoint=mid, low="blue", mid="white",
#                      high="red", space ="Lab" )



# # Install
# install.packages("wesanderson")
# # Load
library(wesanderson)

mp<-mp+scale_color_manual(values=wes_palette(name="Cavalcanti1"))