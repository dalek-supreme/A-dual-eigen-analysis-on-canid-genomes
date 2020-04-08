load('snp.svd.Rdata')

# num_layers <- 10
# for (layer in seq(num_layers)){
#     component <- snp.svd$d[layer] * snp.svd$u[,layer] %*% t(snp.svd$v[,layer]) 
# }

setwd('~/observations')

pdf(file=paste('Histogram-singlar-values.pdf',sep=''))
hist(snp.svd$d)
dev.off()
num_layers <- 10
for (layer in seq(num_layers)){
    pdf(file=paste('Histogram-left-',layer,'.pdf',sep=''))
    hist(snp.svd$u[,layer])
    dev.off()
    pdf(file=paste('Histogram-right-',layer,'.pdf',sep=''))
    hist(snp.svd$v[,layer])
    dev.off()
}


snp.svd.stats <- array(0,c(num_layers,4))
colnames(snp.svd.stats) <- c('left average','left variance','right average','right variance')
rownames(snp.svd.stats) <- seq(num_layers)
for (layer in seq(num_layers)) {
    snp.svd.stats[layer,'left average'] <- mean(snp.svd$u[,layer])
    snp.svd.stats[layer,'left variance'] <- var(snp.svd$u[,layer])
    snp.svd.stats[layer,'right average'] <- mean(snp.svd$v[,layer])
    snp.svd.stats[layer,'right variance'] <- var(snp.svd$v[,layer])
}

snp.svd.stats <- snp.svd.stats[2:num_layers,]

pdf(file='left-average.pdf')
plot(snp.svd.stats[,'left average'])
dev.off()

pdf(file='left-variance.pdf')
plot(snp.svd.stats[,'left variance'])
dev.off()

pdf(file='right-average.pdf')
plot(snp.svd.stats[,'right average'])
dev.off()

pdf(file='right-variance.pdf')
plot(snp.svd.stats[,'right variance'])
dev.off()

