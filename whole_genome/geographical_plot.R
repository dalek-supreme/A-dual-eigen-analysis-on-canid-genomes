library(gdsfmt)
library(SNPRelate)

# for the 127-sample dataset first 
genome <- snpgdsOpen("127_biallelic.gds")
sample.id <- read.gdsn(index.gdsn(genome, "sample.id"))


load('snp.svd.Rdata')
source('DualEigen_utils.R')
## Get main positions:
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

setwd("~/enrichment2/demographic/")
for (layer in seq(num_layers)) {
    demographic.output <- demographic.data[main.pos.x[[layer]],]
    demographic.output$Loading <- snp.svd$u[main.pos.x[[layer]],layer]
    demographic.output <- demographic.output[order((demographic.output$Loading)^2,decreasing=T),]
    write.csv(demographic.output,file=(paste('demographic-',layer,'.csv',sep='')))
}

setwd("~/enrichment2/map")
for (layer in 2:num_layers) {
    demographic.map <- demographic.data[main.pos.x[[layer]],]
    demographic.map <- demographic.map[demographic.map$SampleID!='Lcu2_Pastora',]
    mapworld<-borders("world",colour = "gray50",fill="white")
    mp<-ggplot()+mapworld+ylim(-60,90)
    mp<-mp+geom_point(aes(x=as.numeric(demographic.map$Longitude),
                          y=as.numeric(demographic.map$Latitude)),color='red',size=5,alpha=0.2)
    #mp <- mp+theme(legend.text = element_text('Loading Positivity'))
    pdf(file=paste('map-',layer,'.pdf',sep=''),width = 10,height = 5)
    mp
    dev.off()
}

