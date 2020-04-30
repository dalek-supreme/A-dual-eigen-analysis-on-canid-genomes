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
                        y=as.numeric(cluster.map$Latitude),color=as.factor(cluster.map$ClusterNumber),size=5,alpha=0.2))
# mid<-mean(unique(cluster.map$ClusterNumber))
# mp<-mp+scale_color_gradient2(midpoint=mid, low="blue", mid="white",
#                      high="red", space ="Lab" )



# # Install
# install.packages("wesanderson")
# # Load
library(wesanderson)

mp<-mp+scale_color_manual(values=wes_palette(name="Cavalcanti1"))


library(treeio)
library(ggtree)
snp.tree <- read.tree("snphylo.output.ml.tree")
snp.tree.plot <- ggtree(snp.tree, layout="circular") + geom_tiplab(aes(angle=angle), color='purple',hjust = -.2,size=1.5)+ xlim(0, 0.55) 

setwd('~')
load('snp.svd.Rdata')
source('DualEigen_utils.R')
num_layers <- 10
u.cutoff <- 0.8 -> v.cutoff
main.pos.x <- list()
main.pos.y <- list()
for (layer in seq(num_layers)) {
    u <- filter.loadings(snp.svd$u[,layer],method='2',norm.cutoff=u.cutoff)
    v <- filter.loadings(snp.svd$v[,layer],method='2',norm.cutoff=v.cutoff)
    main.pos.x[[layer]] <- nonzero.pos(u)
    main.pos.y[[layer]] <- nonzero.pos(v)
}

# svd.cluster <- array('black',c(127,10))
# for (layer in 1:10){
#     for (j in seq_along(main.pos.x[[layer]])){
#         svd.cluster[main.pos.x[[layer]][j],layer] <- 'red'
#     }
# }
setwd('tree-dualeigen')
for (layer in 2:10){
svd.cluster2 <- data.frame(
    SampleID = demographic.data$SampleID,
    color = array('grey',length(demographic.data$SampleID))
)

for (j in seq_along(main.pos.x[[layer]])){
        svd.cluster2[main.pos.x[[layer]][j],]$color <- 'red'
}

snp.tree$tip.label[snp.tree$tip.label=='Lcu2_Pasto']<- 'Lcu2_Pastora'
svd.color2 <- array()
for (i in seq_along(svd.cluster2$color)){
    svd.color2[i] <- svd.cluster2[svd.cluster2$SampleID==snp.tree$tip.label[i],]$color
}
pdf(paste('tree-dualeigen-',layer,'.pdf',sep=''))
plot(snp.tree, tip.color=svd.color2, type="fan",cex=0.5,x.lim=0.2)
dev.off()
}
