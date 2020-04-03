num_layers <- 10
demographic.data <- read.csv('demographic.csv')

demographic.index<-0
for (i in seq_along(sample.id)) {
    demographic.index[i]<-which(demographic.data$SampleID==sample.id[i])
}

demographic.data <- demographic.data[demographic.index,]

setwd("enrichment/demographic/")
for (layer in num_layers) {
    canid.order <- order(snp.svd$u[,layer])
    demographic.output <- demographic.data[canid.order,]
    write.csv(demographic.output,file=(paste('demographic_layer_',layer,'.csv',sep='')))
}

