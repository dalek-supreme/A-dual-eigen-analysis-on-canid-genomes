num_layers <- 10
demographic.data <- read.csv('demographic.csv')

demographic.index<-0
for (i in seq_along(sample.id)) {
    demographic.index[i]<-which(demographic.data$SampleID==sample.id[i])
}

demographic.data <- demographic.data[demographic.index,]

library(ggplot2)
library(ggmap)
library(sp)
library(maptools)
library(maps)

setwd("enrichment/demographic/")
for (layer in seq(num_layers)) {
    canid.order <- order(snp.svd$u[,layer])
    demographic.output <- demographic.data[canid.order,]
    write.csv(demographic.output,file=(paste('demographic_layer_',layer,'.csv',sep='')))
}

setwd("enrichment_whole_genome/map")
for (layer in 2:num_layers) {
    demographic.map <- demographic.data
    positivity <- sapply(snp.svd$u[,layer],
                                        function(x) {
                                            if (x<0) return(1)
                                            else if (x>=0) return(2)
                                        })
    demographic.map$color <- factor(positivity,labels=c('negative','positive'))
    demographic.map <- demographic.map[demographic.data$SampleID!='Lcu2_Pastora',]
    mapworld<-borders("world",colour = "gray50",fill="white")
    mp<-ggplot()+mapworld+ylim(-60,90)
    mp<-mp+geom_point(aes(x=as.numeric(demographic.map$Longitude),
                          y=as.numeric(demographic.map$Latitude),color=demographic.map$color))
    #mp <- mp+theme(legend.text = element_text('Loading Positivity'))
    pdf(file=paste('map_layer_',layer,'.pdf',sep=''),width = 10,height = 5)
    mp
    dev.off()
}