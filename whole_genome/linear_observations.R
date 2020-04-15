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

num_obs <- 20
snp.svd.orders <- list()
snp.svd.orders$left <- array(0,c(num_layers,num_obs))
snp.svd.orders$right <- array(0,c(num_layers,num_obs))

for (layer in seq(num_layers)) {
    snp.svd.orders$left[layer,] <- order(snp.svd$u[,layer])[seq(num_obs)]
    snp.svd.orders$right[layer,] <- order(snp.svd$v[,layer])[seq(num_obs)]
}

num_obs <- 20
snp.svd.orders.end <- list()
snp.svd.orders.end$left <- array(0,c(num_layers,num_obs))
snp.svd.orders.end$right <- array(0,c(num_layers,num_obs))

for (layer in seq(num_layers)) {
    snp.svd.orders.end$left[layer,] <- order(snp.svd$u[,layer])[length(snp.svd$u[,layer])+1-seq(num_obs)]
    snp.svd.orders.end$right[layer,] <- order(snp.svd$v[,layer])[length(snp.svd$v[,layer])+1-seq(num_obs)]
}

num_obs <- 20
snp.svd.values <- list()
snp.svd.values$left <- array(0,c(num_layers,num_obs))
snp.svd.values$right <- array(0,c(num_layers,num_obs))

for (layer in seq(num_layers)) {
    snp.svd.values$left[layer,] <- sort(snp.svd$u[,layer])[seq(num_obs)]
    snp.svd.values$right[layer,] <- sort(snp.svd$v[,layer])[seq(num_obs)]
}

num_obs <- 20
snp.svd.values.end <- list()
snp.svd.values.end$left <- array(0,c(num_layers,num_obs))
snp.svd.values.end$right <- array(0,c(num_layers,num_obs))

for (layer in seq(num_layers)) {
    snp.svd.values.end$left[layer,] <- sort(snp.svd$u[,layer])[length(snp.svd$u[,layer])+1-seq(num_obs)]
    snp.svd.values.end$right[layer,] <- sort(snp.svd$v[,layer])[length(snp.svd$v[,layer])+1-seq(num_obs)]
}

cor(order(snp.svd$u[,5])[seq(num_obs)],order(snp.svd$u[,6])[seq(num_obs)])
cor(order(snp.svd$u[,5])[seq(5)],order(snp.svd$u[,6])[seq(5)])
cor(sort(snp.svd$u[,5])[seq(5)],sort(snp.svd$u[,6])[seq(5)])

cor(order(snp.svd$u[,5])[length(snp.svd$u[,5])+1-seq(5)],order(snp.svd$u[,6])[length(snp.svd$u[,6])+1-seq(5)])
cor(sort(snp.svd$u[,5])[length(snp.svd$u[,5])+1-seq(5)],sort(snp.svd$u[,6])[length(snp.svd$u[,6])+1-seq(5)])

layer1_times_1 <- snp.svd$u[,1] %*% t(snp.svd$v[,1])

pdf(file='heatmap1.pdf')
heatmap(layer1_times_1[,1:200])
dev.off()

layer3_times_3 <- sort(snp.svd$u[,3]) %*% t(sort(snp.svd$v[,3]))
pdf(file='heatmap2.pdf')
heatmap(layer3_times_3[,1:200])
dev.off()


layer4_times_4 <- sort(snp.svd$u[,4]) %*% t(sort(snp.svd$v[,4]))

layer4<-array(0,c(127,200))
for (i in 1:(200-1)){
    for (j in 1:127){
        layer4[j,i]<- mean(layer4_times_4[j,((i-1)*1300):(i*1300)])
    }
}
i<-200
for (j in 1:127){
    layer4[j,i]<- mean(layer4_times_4[j,((i-1)*1300):dim(layer4_times_4)[2]])
}


pdf(file='heatmap3.pdf')
heatmap(layer4)
dev.off()

# interesting...

layer5abriv <- array(0,200)
for (i in 1:(200-1)){
    layer5abriv[i]<- mean(sort(snp.svd$v[,5])[((i-1)*1300):(i*1300)])
}
i<-200
layer5abriv[i]<- mean(sort(snp.svd$v[,5])[((i-1)*1300):dim(snp.svd$v)[1]])

layer5_times_5 <- sort(snp.svd$u[,4]) %*% t(sort(layer5abriv))

pdf(file='heatmap4.pdf')
heatmap(layer5_times_5)
dev.off()

a<-rnorm(n=200,sd=100)
b<-rnorm(n=200,sd=1)
pdf(file='heatmap_normal.pdf')
heatmap(a%*%t(b))
dev.off()

layer6ends <- array(0,200)
layer6ends[1:100]<-sort(snp.svd$v[,6])[1:100]
layer6ends[101:200]<-sort(snp.svd$v[,6])[(dim(snp.svd$v)[1]-99):dim(snp.svd$v)[1]]
layer6_times_6 <- sort(snp.svd$u[,6]) %*% t(sort(layer6ends))
pdf(file='heatmap5.pdf')
heatmap(layer6_times_6)
dev.off()


layer1.svd<-svd(snp.svd$u[,1] %*% t(snp.svd$v[,1]))

layer1.sort.svd<-svd(sort(snp.svd$u[,1]) %*% t(sort(snp.svd$v[,1])))

layer1.main.svd<-svd(principle.loadings(snp.svd$u[,1]) %*% t(principle.loadings(snp.svd$v[,1])))

layer2.main.svd<-svd(principle.loadings(snp.svd$u[,2]) %*% t(principle.loadings(snp.svd$v[,2])))

layer3.main.svd<-svd(principle.loadings(snp.svd$u[,2]) %*% t(principle.loadings(snp.svd$v[,2])))

norm2s<-array()
for (i in seq(num_layers)) {
    norm2s[i]<-svd(principle.loadings(snp.svd$u[,i]) %*% t(principle.loadings(snp.svd$v[,i])))$d[1]
    print(norm2s[i])
}

noise.norm2s<-array()
for (i in seq(num_layers)) {
    noise.norm2s[i]<-svd(non.principle.loadings(snp.svd$u[,i]) %*% t(non.principle.loadings(snp.svd$v[,i])))$d[1]
    print(noise.norm2s[i])
}

# compare with norm2s:
norm2s.eval<-array()
for (i in seq(num_layers)) {
    norm2s.eval[i] <- sqrt(sum((snp.svd$u[,i])^2)*sum((snp.svd$v[,i])^2))
    print(norm2s.eval)
}

# validify the result on the 2-norm of a rank-1 matrix
a<-c(1,2,3,4,5)
b<-c(2,24,43,1)
a%*%t(b)
norm(a%*%t(b),type='2')
svd(a%*%t(b))$d[1]

sqrt(sum(a^2)*sum(b^2))

# the operator norm of rank 1 matrix = the product of the norms of the two vectors

noise.norm2s.eval<-array()
for (i in seq(num_layers)) {
    noise.norm2s.eval[i] <- sqrt(sum(non.principle.loadings(snp.svd$u[,i])^2)*sum(non.principle.loadings(snp.svd$v[,i])^2))
    print(noise.norm2s.eval)
}

#  [1] 0.81946857 0.04639882 0.07888633 0.10814524 0.07286147 0.07845373
#  [7] 0.02964882 0.02618803 0.08435135 0.01949680

