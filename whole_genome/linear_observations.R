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

#  [1] 0.81946857 0.04639882 0.07888633 0.10814524 0.07286147 0.07845373
#  [7] 0.02964882 0.02618803 0.08435135 0.01949680



############### pertubation and noise and stability
pdf('plot_pertub.pdf')
plot(diff(snp.svd$d)/snp.svd$d[2:length(snp.svd$d)][2:(length(snp.svd$d))])
dev.off()

diff(snp.svd$d)/snp.svd$d[2:length(snp.svd$d)]->pertub

(snp.svd$d[1:(length(snp.svd$d)-1)]-snp.svd$d[2:length(snp.svd$d)])

# get principle loadings, combine and sum to get approximation matrix, calculate error. (2-norm, Frobenius norm)

a<-array(0,c(10,10))
a[1,3]<-1
a[3,5]<-1
a[4,9]<-1
a[3,4]<-1

b<-svd(a)

b$d[1]*b$u[,1]%*%t(b$v[,1])+b$d[2]*b$u[,2]%*%t(b$v[,2])+b$d[3]*b$u[,3]%*%t(b$v[,3])


# check sparsity
count.0 <- 0
count.1 <- 0
count.2 <- 0

for (i in seq(dim(snp.genotype)[1])){
    for (j in seq(dim(snp.genotype)[2])){
        if (snp.genotype[i,j]==0)
            count.0 <- count.0+1
        else if (snp.genotype[i,j]==1)
            count.1 <- count.1+1
        else if (snp.genotype[i,j]==2)
            count.2 <- count.2+1
        else print(snp.genotype[i,j])
    }
}
print(count.0)
print(count.1)
print(count.2)

count.0+count.1+count.2
dim(snp.genotype)[1]*dim(snp.genotype)[2]

trivial.cols <- 0
for (j in seq(dim(snp.genotype)[2])){
    if (  (!(1 %in% snp.genotype[,j]) & !(0 %in% snp.genotype[,j]))
        | (!(2 %in% snp.genotype[,j]) & !(0 %in% snp.genotype[,j]))
        | (!(1 %in% snp.genotype[,j]) & !(2 %in% snp.genotype[,j]))
        ) {
        trivial.cols <- trivial.cols +1
        print(trivial.cols)
    }
}

trivial.rows <- 0
for (i in seq(dim(snp.genotype)[1])){
    if (  (!(1 %in% snp.genotype[i,]) & !(0 %in% snp.genotype[i,]))
        | (!(2 %in% snp.genotype[i,]) & !(0 %in% snp.genotype[i,]))
        | (!(1 %in% snp.genotype[i,]) & !(2 %in% snp.genotype[i,]))
        ) {
        trivial.rows <- trivial.rows +1
        print(trivial.rows)
    }
}
# no trivial columns or rows
# indicates very well distributed sparsity
# full rank can be a problem...

# is full rank since it has 127 singular values.\

# change coding to 2,1,0

snp.genotype.recoding <- snp.genotype

for (i in seq(dim(snp.genotype.recoding)[1])){
    for (j in seq(dim(snp.genotype.recoding)[2])){
        if (snp.genotype.recoding[i,j]==0)
            snp.genotype.recoding[i,j]<-2
        else if (snp.genotype.recoding[i,j]==2)
            snp.genotype.recoding[i,j]<-0
    }
}


# check sparsity
count.0 <- 0
count.1 <- 0
count.2 <- 0

for (i in seq(dim(snp.genotype.recoding)[1])){
    for (j in seq(dim(snp.genotype.recoding)[2])){
        if (snp.genotype.recoding[i,j]==0)
            count.0 <- count.0+1
        else if (snp.genotype.recoding[i,j]==1)
            count.1 <- count.1+1
        else if (snp.genotype.recoding[i,j]==2)
            count.2 <- count.2+1
        else print(snp.genotype.recoding[i,j])
    }
}
print(count.0)
print(count.1)
print(count.2)


load('snp.recoding.svd.Rdata')
svd.comp <- array(0,c(2,127))
svd.comp[1,] <- snp.svd$d
svd.comp[2,] <- snp.recoding.svd$d

sqrt(sum((snp.svd$u[,1]-snp.recoding.svd$u[,1])^2))
sqrt(sum((snp.svd$v[,1]-snp.recoding.svd$v[,1])^2))
sqrt(sum((snp.svd$v[,1])^2))
sqrt(sum((snp.recoding.svd$v[,1])^2))

sum(snp.svd$v[,1]*snp.recoding.svd$v[,1])
acos(sum(snp.svd$v[,1]*snp.recoding.svd$v[,1]))

acos(sum(snp.svd$v[,1]*snp.recoding.svd$v[,1]))/pi*180
acos(sum(snp.svd$v[,2]*snp.recoding.svd$v[,2]))/pi*180
acos(sum(snp.svd$v[,3]*snp.recoding.svd$v[,3]))/pi*180

for (i in seq(10)) {
    print(acos(sum(snp.svd$v[,i]*snp.recoding.svd$v[,i]))/pi*180)
}
for (i in seq(10)) {
    print(acos(sum(snp.svd$u[,i]*snp.recoding.svd$u[,i]))/pi*180)
}

for (i in seq(10)){
    print(sum(order(snp.svd$u[,i])-order(snp.recoding.svd$u[,i])))
}
for (i in seq(10)){
    print(sum(order(snp.svd$v[,i])-order(snp.recoding.svd$v[,i])))
}

# interestingly, no order has changed.

A1<-snp.recoding.svd$d[1]*snp.recoding.svd$u[,1]%*%t(snp.recoding.svd$v[,1])

A1.rowmax<-array()
for (i in seq(dim(A1)[1])){
    A1.rowmax[i]<-max(A1[i,])
}

A1.colmax<-array()
for (j in seq(dim(A1)[2])){
    A1.colmax[j]<-max(A1[,j])
}


snp.genotype.recoding.rowmax<-array()
for (i in seq(dim(snp.genotype.recoding)[1])){
    snp.genotype.recoding.rowmax[i]<-max(snp.genotype.recoding[i,])
}

snp.genotype.recoding.colmax<-array()
for (j in seq(dim(snp.genotype.recoding)[2])){
    snp.genotype.recoding.colmax[j]<-max(snp.genotype.recoding[,j])
}


setwd('~/observations/recoding')

pdf(file=paste('Histogram-singlar-values.pdf',sep=''))
hist(snp.recoding.svd$d)
dev.off()
num_layers <- 10
for (layer in seq(num_layers)){
    pdf(file=paste('Histogram-left-',layer,'.pdf',sep=''))
    hist(snp.recoding.svd$u[,layer])
    dev.off()
    pdf(file=paste('Histogram-right-',layer,'.pdf',sep=''))
    hist(snp.recoding.svd$v[,layer])
    dev.off()
}


A1<-snp.recoding.svd$d[1]*snp.recoding.svd$u[,1]%*%t(snp.recoding.svd$v[,1])
hist(A1)
A2<-snp.recoding.svd$d[2]*snp.recoding.svd$u[,2]%*%t(snp.recoding.svd$v[,2])
hist(A2)

sort(abs(A1))[(length(A1)-100):length(A1)]
sort(abs(A2))[(length(A2)-100):length(A2)]

v2 <- filter.loadings(snp.svd$v[,2],method='2',norm.cutoff=0.9)
count <- 0
for (val in v2){
    if (val != 0) {
        count <- count +1
    }
}
count
count/length(v2)

for (layer in seq(num_layers)){
    v2 <- filter.loadings(snp.svd$v[,layer],method='2',norm.cutoff=0.9)
    count <- 0
    for (val in v2){
        if (val != 0) {
            count <- count +1
        }
    }
    print(count)
    print(count/length(v2))
}

for (layer in seq(num_layers)){
    v2 <- filter.loadings(snp.svd$u[,layer],method='2',norm.cutoff=0.9)
    count <- 0
    for (val in v2){
        if (val != 0) {
            count <- count +1
        }
    }
    print(count)
    print(count/length(v2))
}






num_layers<-10
B <- array(0,dim(snp.genotype))
for (layer in seq(num_layers)){
    u <- filter.loadings(snp.svd$u[,layer],method='2',norm.cutoff=0.99)
    v <- filter.loadings(snp.svd$v[,layer],method='2',norm.cutoff=0.99)
    B <- B + snp.svd$d[layer]*u %*% t(v)
}


err <- norm(snp.genotype-B, type='2')

A1 <- snp.svd$u[,1]%*%t(snp.svd$v[,1])
u <- filter.loadings(snp.svd$u[,1],method='2',norm.cutoff=0.9)
v <- filter.loadings(snp.svd$v[,1],method='2',norm.cutoff=0.9)
B1 <- u %*% t(v)

norm(A1-B1,type='2')

A2 <- snp.svd$u[,2]%*%t(snp.svd$v[,2])
u <- filter.loadings(snp.svd$u[,2],method='2',norm.cutoff=0.9)
v <- filter.loadings(snp.svd$v[,2],method='2',norm.cutoff=0.9)
B2 <- u %*% t(v)

norm(A2-B2,type='2')

A3 <- snp.svd$u[,3]%*%t(snp.svd$v[,3])
u <- filter.loadings(snp.svd$u[,3],method='2',norm.cutoff=0.9)
v <- filter.loadings(snp.svd$v[,3],method='2',norm.cutoff=0.9)
B3 <- u %*% t(v)

norm(A3-B3,type='2')

A4 <- snp.svd$u[,4]%*%t(snp.svd$v[,4])
u <- filter.loadings(snp.svd$u[,4],method='2',norm.cutoff=0.9)
v <- filter.loadings(snp.svd$v[,4],method='2',norm.cutoff=0.9)
B4 <- u %*% t(v)

norm(A4-B4,type='2')


for (layer in seq(num_layers)){
    u <- filter.loadings(snp.svd$u[,layer],method='2',norm.cutoff=0.7)
    v <- filter.loadings(snp.svd$v[,layer],method='2',norm.cutoff=0.7)
    count0 <- 0
    for (i in u) {
        if (i!=0) count0 <- count0+1
    }
    print(paste(layer,'u',count0))
    count0 <- 0
    for (i in v) {
        if (i!=0) count0 <- count0+1
    }
    print(paste(layer,'v',count0))
}

num_layers<-10
B <- snp.svd$d[1]*snp.svd$u[,1]%*%t(snp.svd$v[,1])
for (layer in 2:num_layers){
    u <- filter.loadings(snp.svd$u[,layer],method='2',norm.cutoff=0.7)
    v <- filter.loadings(snp.svd$v[,layer],method='2',norm.cutoff=0.7)
    B <- B + snp.svd$d[layer]*u %*% t(v)
}
err <- norm(snp.genotype-B, type='2')

for (layer in seq(num_layers)){
    u <- filter.loadings(snp.svd$u[,layer],method='2',norm.cutoff=0.6)
    v <- filter.loadings(snp.svd$v[,layer],method='2',norm.cutoff=0.6)
    count0 <- 0
    for (i in u) {
        if (i!=0) count0 <- count0+1
    }
    print(paste(layer,'u',count0))
    count0 <- 0
    for (i in v) {
        if (i!=0) count0 <- count0+1
    }
    print(paste(layer,'v',count0))
}

num_layers<-10
B <- snp.svd$d[1]*snp.svd$u[,1]%*%t(snp.svd$v[,1])
for (layer in 2:num_layers){
    u <- filter.loadings(snp.svd$u[,layer],method='2',norm.cutoff=0.6)
    v <- filter.loadings(snp.svd$v[,layer],method='2',norm.cutoff=0.6)
    B <- B + snp.svd$d[layer]*u %*% t(v)
}
err <- norm(snp.genotype-B, type='2')


for (layer in seq(num_layers)){
    u <- filter.loadings(snp.svd$u[,layer],method='2',norm.cutoff=0.5)
    v <- filter.loadings(snp.svd$v[,layer],method='2',norm.cutoff=0.5)
    count0 <- 0
    for (i in u) {
        if (i!=0) count0 <- count0+1
    }
    print(paste(layer,'u',count0))
    count0 <- 0
    for (i in v) {
        if (i!=0) count0 <- count0+1
    }
    print(paste(layer,'v',count0))
}

num_layers<-10
B <- snp.svd$d[1]*snp.svd$u[,1]%*%t(snp.svd$v[,1])
for (layer in 2:num_layers){
    u <- filter.loadings(snp.svd$u[,layer],method='2',norm.cutoff=0.5)
    v <- filter.loadings(snp.svd$v[,layer],method='2',norm.cutoff=0.5)
    B <- B + snp.svd$d[layer]*u %*% t(v)
}
err <- norm(snp.genotype-B, type='2')

cutoff <- 0.4

num_layers<- 8
u.cutoff <- 0.8
v.cutoff <- 0.8
for (layer in seq(num_layers)){
    u <- filter.loadings(snp.svd$u[,layer],method='2',norm.cutoff=u.cutoff)
    v <- filter.loadings(snp.svd$v[,layer],method='2',norm.cutoff=v.cutoff)
    count0 <- 0
    for (i in u) {
        if (i!=0) count0 <- count0+1
    }
    print(paste(layer,'u',count0))
    count0 <- 0
    for (i in v) {
        if (i!=0) count0 <- count0+1
    }
    print(paste(layer,'v',count0))
}

B <- snp.svd$d[1]*snp.svd$u[,1]%*%t(snp.svd$v[,1])
for (layer in 2:num_layers){
    u <- filter.loadings(snp.svd$u[,layer],method='2',norm.cutoff=u.cutoff)
    v <- filter.loadings(snp.svd$v[,layer],method='2',norm.cutoff=v.cutoff)
    B <- B + snp.svd$d[layer]*u %*% t(v)
}
err <- norm(snp.genotype-B, type='2')


for (i in 1:10) {
    u.cutoff <- 1-0.02*i -> v.cutoff
    print(u.cutoff)
    for (layer in seq(num_layers)){
        u <- filter.loadings(snp.svd$u[,layer],method='2',norm.cutoff=u.cutoff)
        v <- filter.loadings(snp.svd$v[,layer],method='2',norm.cutoff=v.cutoff)
        count0 <- 0
        for (i in u) {
            if (i!=0) count0 <- count0+1
        }
        print(paste(layer,'u',count0))
        count0 <- 0
        for (i in v) {
            if (i!=0) count0 <- count0+1
        }
        print(paste(layer,'v',count0))
    }

    num_layers<-10
    B <- snp.svd$d[1]*snp.svd$u[,1]%*%t(snp.svd$v[,1])
    for (layer in 2:num_layers){
        u <- filter.loadings(snp.svd$u[,layer],method='2',norm.cutoff=u.cutoff)
        v <- filter.loadings(snp.svd$v[,layer],method='2',norm.cutoff=v.cutoff)
        B <- B + snp.svd$d[layer]*u %*% t(v)
    }
    err <- norm(snp.genotype-B, type='2')
    print(err)
}

u.cutoff <- 0.8 -> v.cutoff
for (num_layers in 1:10){
    B <- snp.svd$d[1]*snp.svd$u[,1]%*%t(snp.svd$v[,1])
    for (layer in 2:num_layers){
        u <- filter.loadings(snp.svd$u[,layer],method='2',norm.cutoff=u.cutoff)
        v <- filter.loadings(snp.svd$v[,layer],method='2',norm.cutoff=v.cutoff)
        B <- B + snp.svd$d[layer]*u %*% t(v)
    }
    err <- norm(snp.genotype-B, type='2')
    print(paste(num_layers,err))
}

C <- array(2,dim(snp.genotype))
norm(C,type='2')
2*sqrt(dim(snp.genotype)[1]*dim(snp.genotype)[2])

norm(C-snp.svd$d[1]*snp.svd$u[,1]%*%t(snp.svd$v[,1]),type='2')
# 1604.177 

norm(C-snp.svd$d[1]*snp.svd$u[,1]%*%t(snp.svd$v[,1])
      -snp.svd$d[2]*snp.svd$u[,2]%*%t(snp.svd$v[,2])
      ,type='2')
# 1604.431

norm(C-snp.svd$d[1]*snp.svd$u[,1]%*%t(snp.svd$v[,1])
      -snp.svd$d[2]*snp.svd$u[,2]%*%t(snp.svd$v[,2])
      -snp.svd$d[3]*snp.svd$u[,3]%*%t(snp.svd$v[,3])
      ,type='2')
# 1604.435

order(snp.old.svd$u[,1])==order(snp.recoding.svd$u[,1])

num_layers<- 8
u.cutoff <- 0.8
v.cutoff <- 0.8
for (layer in seq(num_layers)){
    u <- filter.loadings(snp.svd$u[,layer],method='2',norm.cutoff=u.cutoff)
    v <- filter.loadings(snp.svd$v[,layer],method='2',norm.cutoff=v.cutoff)
    print(paste(layer,'u'))
    print(nonzero.pos(u))
}
for (layer in seq(num_layers)){
    u <- filter.loadings(snp.svd$u[,layer],method='2',norm.cutoff=u.cutoff)
    v <- filter.loadings(snp.svd$v[,layer],method='2',norm.cutoff=v.cutoff)
    print(paste(layer,'v'))
    print(nonzero.pos(v))
}

cor_mat <- array(0,c(num_layers,num_layers))
for (i in seq(num_layers-1)){
    for (j in (i+1):num_layers){
        u1 <- filter.loadings(snp.svd$u[,i],method='2',norm.cutoff=u.cutoff)
        u2 <- filter.loadings(snp.svd$u[,j],method='2',norm.cutoff=u.cutoff)
        cor_mat[i,j]<-cor(u1,u2)
    }
}

cor_mat <- array(0,c(num_layers,num_layers))
for (i in seq(num_layers-1)){
    for (j in (i+1):num_layers){
        u1 <- filter.loadings(snp.svd$u[,i],method='2',norm.cutoff=u.cutoff)
        u2 <- filter.loadings(snp.svd$u[,j],method='2',norm.cutoff=u.cutoff)
        cor_mat[i,j]<-cor(order(u1),order(u2))
    }
}

inter_mat <- array(0,c(num_layers,num_layers))
for (i in seq(num_layers-1)){
    for (j in (i+1):num_layers){
        u1 <- filter.loadings(snp.svd$u[,i],method='2',norm.cutoff=u.cutoff)
        u2 <- filter.loadings(snp.svd$u[,j],method='2',norm.cutoff=u.cutoff)
        inter_mat[i,j]<-length(intersect(nonzero.pos(u1),nonzero.pos(u2)))
    }
}

u.inter_rate_mat <- array(0,c(num_layers,num_layers))
for (i in seq(num_layers-1)){
    for (j in (i+1):num_layers){
        u1 <- filter.loadings(snp.svd$u[,i],method='2',norm.cutoff=u.cutoff)
        u2 <- filter.loadings(snp.svd$u[,j],method='2',norm.cutoff=u.cutoff)
        u.inter_rate_mat[i,j]<-length(intersect(nonzero.pos(u1),nonzero.pos(u2)))/length(union(nonzero.pos(u1),nonzero.pos(u2)))
    }
}
v.inter_rate_mat <- array(0,c(num_layers,num_layers))
for (i in seq(num_layers-1)){
    for (j in (i+1):num_layers){
        v1 <- filter.loadings(snp.svd$v[,i],method='2',norm.cutoff=v.cutoff)
        v2 <- filter.loadings(snp.svd$v[,j],method='2',norm.cutoff=v.cutoff)
        v.inter_rate_mat[i,j]<-length(intersect(nonzero.pos(v1),nonzero.pos(v2)))/length(union(nonzero.pos(v1),nonzero.pos(v2)))
    }
}




# MLE of bernouli probability:
count.0 <- 0
count.1 <- 0
count.2 <- 0

for (i in seq(dim(snp.genotype)[1])){
    for (j in seq(dim(snp.genotype)[2])){
        if (snp.genotype[i,j]==0)
            count.0 <- count.0+1
        else if (snp.genotype[i,j]==1)
            count.1 <- count.1+1
        else if (snp.genotype[i,j]==2)
            count.2 <- count.2+1
        else print(snp.genotype[i,j])
    }
}

p <- (count.1+2*count.0)/(2*(count.0+count.1+count.2))
n <- 2*(count.0+count.1+count.2)
derivative5 <- snp.svd$u[,5]%*%t(snp.svd$v[,5])

k <- (snp.svd$d[5]-snp.svd$d[6])/(mean(derivative5)*(2*count.0+count.1)/(count.0+count.1))
2*(1-pnorm(k/sqrt(2*n*p*(1-p))))
# =0

k <- (snp.svd$d[5]-snp.svd$d[6])/(max(derivative5)*(2*count.0+count.1)/(count.0+count.1))
2*(1-pnorm(k/sqrt(2*n*p*(1-p))))

plot(sort(derivative5,decreasing=T)[1:10000])


num_layers<- 8
u.thrshld <- 3/sqrt(length(snp.svd$u[,1]))
v.thrshld <- 3/sqrt(length(snp.svd$v[,1]))
for (layer in seq(num_layers)){
    u <- filter.loadings(snp.svd$u[,layer],method='threshold',threshold=u.thrshld)
    v <- filter.loadings(snp.svd$v[,layer],method='threshold',threshold=v.thrshld)
    print(paste(layer,'u'))
    print(nonzero.pos(u))
}
for (layer in seq(num_layers)){
    u <- filter.loadings(snp.svd$u[,layer],method='threshold',threshold=u.thrshld)
    v <- filter.loadings(snp.svd$v[,layer],method='threshold',threshold=v.thrshld)
    print(paste(layer,'v'))
    print(nonzero.pos(v))
}


num_layers<- 8
u.cutoff <- 0.8
v.cutoff <- 0.8
for (layer in seq(num_layers)){
    u <- filter.loadings(snp.svd$u[,layer],method='2',norm.cutoff=u.cutoff)
    #v <- filter.loadings(snp.svd$v[,layer],method='2',norm.cutoff=v.cutoff)
    print(paste(layer,'u'))
    print(u[nonzero.pos(u)])
}
for (layer in seq(num_layers)){
    #u <- filter.loadings(snp.svd$u[,layer],method='2',norm.cutoff=u.cutoff)
    v <- filter.loadings(snp.svd$v[,layer],method='2',norm.cutoff=v.cutoff)
    print(paste(layer,'v'))
    print(v[nonzero.pos(v)])
}

mutation.count <- array(0,dim(snp.genotype)[1])
for (i in seq(dim(snp.genotype)[1])){
    mutation.count[i] <- sum(snp.genotype[i,])
}

mutation.count <- array(0,c(3,dim(snp.genotype)[1]))
rownames(mutation.count) <- c(0,1,2)
for (i in seq(dim(snp.genotype)[1])){
    for (j in seq(dim(snp.genotype)[2])) {
        if (snp.genotype[i,j]==0)
            mutation.count[1,] <- mutation.count[1,i]+1
        else if (snp.genotype[i,j]==1)
            mutation.count[2,i] <- mutation.count[2,i]+1
        else if (snp.genotype[i,j]==2)
            mutation.count[3,i] <- mutation.count[3,i]+1
        else print(snp.genotype[i,j])
    }
}

# see if the 'main positions' intersect
u.cutoff <- 0.8 -> v.cutoff
num_layers <- 10
main.pos.x <- list()
main.pos.y <- list()
for (layer in seq(num_layers)) {
    u <- filter.loadings(snp.svd$u[,layer],method='2',norm.cutoff=u.cutoff)
    v <- filter.loadings(snp.svd$v[,layer],method='2',norm.cutoff=v.cutoff)
    main.pos.x[[layer]] <- nonzero.pos(u)
    main.pos.y[[layer]] <- nonzero.pos(v)
}

plot()
