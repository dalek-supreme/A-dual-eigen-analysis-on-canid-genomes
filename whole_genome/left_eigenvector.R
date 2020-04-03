demographic.data <- read.csv('demographic.csv')

demographic.index<-0
for (i in seq_along(sample.id)) {
    demographic.index[i]<-which(demographic.data$SampleID==sample.id[i])
}

demographic.data <- demographic.data[demographic.index,]

for (layer in num_layers) {
    canid.order <- order(snp.svd$u[,layer])
    demographic.data[canid.order]
}