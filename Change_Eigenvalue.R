load('snp.genotype.Rdata')
load('snp.svd.Rdata')
source('DualEigen_utils.R')
u.cutoff <- 0.8 -> v.cutoff
num_layers <- 10
B <- snp.svd$d[1]*snp.svd$u[,1]%*%t(snp.svd$v[,1])
for (layer in 2:num_layers){
    u <- filter.loadings(snp.svd$u[,layer],method='2',norm.cutoff=u.cutoff)
    v <- filter.loadings(snp.svd$v[,layer],method='2',norm.cutoff=v.cutoff)
    B <- B + snp.svd$d[layer]*u %*% t(v)
}
err <- norm(snp.genotype-B, type='2')

length.pos.x <- array(0,num_layers)
length.pos.y <- array(0,num_layers)
for (layer in seq(num_layers)) {
    u <- filter.loadings(snp.svd$u[,layer],method='2',norm.cutoff=u.cutoff)
    v <- filter.loadings(snp.svd$v[,layer],method='2',norm.cutoff=v.cutoff)
    length.pos.x[layer] <- length(nonzero.pos(u))
    length.pos.y[layer] <- length(nonzero.pos(v))
}


# distribution within the entire matrix
count.all <- array(0,3)
names(count.all) <- c(0,1,2)
for (i in seq(dim(snp.genotype)[1])){
    for (j in seq(dim(snp.genotype)[2])){
        if (snp.genotype[i,j]==0)
            count.all[1] <- count.all[1]+1
        else if (snp.genotype[i,j]==1)
            count.all[2] <- count.all[2]+1
        else if (snp.genotype[i,j]==2)
            count.all[3] <- count.all[3]+1
        else print(snp.genotype[i,j])
    }
}
probab.all <- count.all/sum(count.all)

# distributions within the main area
num_layers <- 10
probab <- array(0,c(num_layers,3))
for (layer in seq(num_layers)) {
    u <- filter.loadings(snp.svd$u[,layer],method='2',norm.cutoff=u.cutoff)
    v <- filter.loadings(snp.svd$v[,layer],method='2',norm.cutoff=v.cutoff)
    main.pos.x <- nonzero.pos(u)
    main.pos.y <- nonzero.pos(v)
    main.matrix <- snp.genotype[main.pos.x,main.pos.y]

    count <- array(0,3)
    names(count) <- c(0,1,2)
    for (i in seq(dim(main.matrix)[1])){
        for (j in seq(dim(main.matrix)[2])){
            if (main.matrix[i,j]==0)
                count[1] <- count[1]+1
            else if (main.matrix[i,j]==1)
                count[2] <- count[2]+1
            else if (main.matrix[i,j]==2)
                count[3] <- count[3]+1
            else print(main.matrix[i,j])
        }
    }
    probab[layer,] <- count/sum(count)
}

# placebo: random selection
random.pos.x <- random.pos(dim(snp.genotype)[1],length(main.pos.x))
random.pos.y <- random.pos(dim(snp.genotype)[2],length(main.pos.y))
random.matrix <- snp.genotype[random.pos.x,random.pos.y]
count.random <- array(0,3)
names(count.random) <- c(0,1,2)
for (i in seq(dim(random.matrix)[1])){
    for (j in seq(dim(random.matrix)[2])){
        if (random.matrix[i,j]==0)
            count.random[1] <- count.random[1]+1
        else if (random.matrix[i,j]==1)
            count.random[2] <- count.random[2]+1
        else if (random.matrix[i,j]==2)
            count.random[3] <- count.random[3]+1
        else print(random.matrix[i,j])
    }
}
probab.random <- count.random/sum(count.random)
# significantly different


### for layer 2:
layer <- 2
u <- filter.loadings(snp.svd$u[,layer],method='2',norm.cutoff=u.cutoff)
v <- filter.loadings(snp.svd$v[,layer],method='2',norm.cutoff=v.cutoff)
main.pos.x <- nonzero.pos(u)
main.pos.y <- nonzero.pos(v)
main.matrix <- snp.genotype[main.pos.x,main.pos.y]

count <- array(0,3)
names(count) <- c(0,1,2)
for (i in seq(dim(main.matrix)[1])){
    for (j in seq(dim(main.matrix)[2])){
        if (main.matrix[i,j]==0)
            count[1] <- count[1]+1
        else if (main.matrix[i,j]==1)
            count[2] <- count[2]+1
        else if (main.matrix[i,j]==2)
            count[3] <- count[3]+1
        else print(main.matrix[i,j])
    }
}
probab <- count/sum(count)
## replace with random distribution with local distribution:
snp.genotype.alter.2.local <- snp.genotype
#generate
rand.value <- array(runif(dim(main.matrix)[1]*dim(main.matrix)[2]),dim(main.matrix))
rand.genotype <- array(-1,dim(main.matrix))
for (i in seq(dim(main.matrix)[1])){
    for (j in seq(dim(main.matrix)[2])){
        if (rand.value[i,j] < probab[1]) rand.genotype[i,j] <- 0
        else if (rand.value[i,j] < probab[1]+probab[2]) rand.genotype[i,j] <- 1
        else rand.genotype[i,j] <- 2
    }
}
#check
rand.count <- array(0,3)
names(rand.count) <- c(0,1,2)
for (i in seq(dim(rand.genotype)[1])){
    for (j in seq(dim(rand.genotype)[2])){
        if (rand.genotype[i,j]==0)
            rand.count[1] <- rand.count[1]+1
        else if (rand.genotype[i,j]==1)
            rand.count[2] <- rand.count[2]+1
        else if (rand.genotype[i,j]==2)
            rand.count[3] <- rand.count[3]+1
        else print(rand.genotype[i,j])
    }
}
rand.probab <- rand.count/sum(rand.count)
(rand.probab - probab)/probab

# replace
snp.genotype.alter.2.local[main.pos.x,main.pos.y] <- rand.genotype
#svd
snp.svd.alter.2.local <- svd(snp.genotype.alter.2.local)
#comparison
(snp.svd.alter.2.local$d-snp.svd$d)/snp.svd$d
#plot
pdf('snp.svd.alter.2.local.pdf')
plot(((snp.svd.alter.2.local$d-snp.svd$d)/snp.svd$d)[1:num_layers])
dev.off()

## replace with random distribution with global distribution:
snp.genotype.alter.2.global <- snp.genotype
#generate
rand.value <- array(runif(dim(main.matrix)[1]*dim(main.matrix)[2]),dim(main.matrix))
rand.genotype <- array(-1,dim(main.matrix))
for (i in seq(dim(main.matrix)[1])){
    for (j in seq(dim(main.matrix)[2])){
        if (rand.value[i,j] < probab.all[1]) rand.genotype[i,j] <- 0
        else if (rand.value[i,j] < probab.all[1]+probab.all[2]) rand.genotype[i,j] <- 1
        else rand.genotype[i,j] <- 2
    }
}
#check
rand.count <- array(0,3)
names(rand.count) <- c(0,1,2)
for (i in seq(dim(rand.genotype)[1])){
    for (j in seq(dim(rand.genotype)[2])){
        if (rand.genotype[i,j]==0)
            rand.count[1] <- rand.count[1]+1
        else if (rand.genotype[i,j]==1)
            rand.count[2] <- rand.count[2]+1
        else if (rand.genotype[i,j]==2)
            rand.count[3] <- rand.count[3]+1
        else print(rand.genotype[i,j])
    }
}
rand.probab <- rand.count/sum(rand.count)
(rand.probab - probab.all)/probab.all

# replace
snp.genotype.alter.2.global[main.pos.x,main.pos.y] <- rand.genotype
#svd
snp.svd.alter.2.global <- svd(snp.genotype.alter.2.global)
#comparison
(snp.svd.alter.2.global$d-snp.svd$d)/snp.svd$d
#plot
pdf('snp.svd.alter.2.global.pdf')
plot(((snp.svd.alter.2.global$d-snp.svd$d)/snp.svd$d)[1:num_layers])
dev.off()

# ## add 1
snp.genotype.alter.2.add1 <- snp.genotype
snp.genotype.alter.2.add1[main.pos.x,main.pos.y] <- snp.genotype.alter.2.add1[main.pos.x,main.pos.y]+1
#svd
snp.svd.alter.2.add1 <- svd(snp.genotype.alter.2.add1)
#comparison
(snp.svd.alter.2.add1$d-snp.svd$d)/snp.svd$d
#plot
pdf('snp.svd.alter.2.add1.pdf')
plot(((snp.svd.alter.2.add1$d-snp.svd$d)/snp.svd$d)[1:num_layers])
dev.off()

## placebo
# replace elsewhere a same size matrix
# get a random position
random.pos.x <- random.pos(dim(snp.genotype)[1],length(main.pos.x))
random.pos.y <- random.pos(dim(snp.genotype)[2],length(main.pos.y))

# replace with a matrix generated by local distribution
#generate
rand.value <- array(runif(dim(main.matrix)[1]*dim(main.matrix)[2]),dim(main.matrix))
rand.genotype <- array(-1,dim(main.matrix))
for (i in seq(dim(main.matrix)[1])){
    for (j in seq(dim(main.matrix)[2])){
        if (rand.value[i,j] < probab[1]) rand.genotype[i,j] <- 0
        else if (rand.value[i,j] < probab[1]+probab[2]) rand.genotype[i,j] <- 1
        else rand.genotype[i,j] <- 2
    }
}
#check
rand.count <- array(0,3)
names(rand.count) <- c(0,1,2)
for (i in seq(dim(rand.genotype)[1])){
    for (j in seq(dim(rand.genotype)[2])){
        if (rand.genotype[i,j]==0)
            rand.count[1] <- rand.count[1]+1
        else if (rand.genotype[i,j]==1)
            rand.count[2] <- rand.count[2]+1
        else if (rand.genotype[i,j]==2)
            rand.count[3] <- rand.count[3]+1
        else print(rand.genotype[i,j])
    }
}
rand.probab <- rand.count/sum(rand.count)
(rand.probab - probab)/probab

# replace
snp.genotype.alter.2.random.local <- snp.genotype
snp.genotype.alter.2.random.local[random.pos.x,random.pos.y] <- rand.genotype

#svd
snp.svd.alter.2.random.local <- svd(snp.genotype.alter.2.random.local)
#comparison
(snp.svd.alter.2.random.local$d-snp.svd$d)/snp.svd$d
#plot
pdf('snp.svd.alter.2.random.local.pdf')
plot(((snp.svd.alter.2.random.local$d-snp.svd$d)/snp.svd$d)[1:num_layers])
dev.off()

# replace with a matrix generated by global distribution
#generate
rand.value <- array(runif(dim(main.matrix)[1]*dim(main.matrix)[2]),dim(main.matrix))
rand.genotype <- array(-1,dim(main.matrix))
for (i in seq(dim(main.matrix)[1])){
    for (j in seq(dim(main.matrix)[2])){
        if (rand.value[i,j] < probab.all[1]) rand.genotype[i,j] <- 0
        else if (rand.value[i,j] < probab.all[1]+probab.all[2]) rand.genotype[i,j] <- 1
        else rand.genotype[i,j] <- 2
    }
}
#check
rand.count <- array(0,3)
names(rand.count) <- c(0,1,2)
for (i in seq(dim(rand.genotype)[1])){
    for (j in seq(dim(rand.genotype)[2])){
        if (rand.genotype[i,j]==0)
            rand.count[1] <- rand.count[1]+1
        else if (rand.genotype[i,j]==1)
            rand.count[2] <- rand.count[2]+1
        else if (rand.genotype[i,j]==2)
            rand.count[3] <- rand.count[3]+1
        else print(rand.genotype[i,j])
    }
}
rand.probab <- rand.count/sum(rand.count)
(rand.probab - probab.all)/probab

# replace
snp.genotype.alter.2.random.global <- snp.genotype
snp.genotype.alter.2.random.global[random.pos.x,random.pos.y] <- rand.genotype

#svd
snp.svd.alter.2.random.global <- svd(snp.genotype.alter.2.random.global)
#comparison
(snp.svd.alter.2.random.global$d-snp.svd$d)/snp.svd$d
#plot
pdf('snp.svd.alter.2.random.global.pdf')
plot(((snp.svd.alter.2.random.global$d-snp.svd$d)/snp.svd$d)[1:num_layers])
dev.off()

# ## add 1
snp.genotype.alter.2.add1.randpos <- snp.genotype
snp.genotype.alter.2.add1.randpos[random.pos.x,random.pos.y] <- snp.genotype.alter.2.add1.randpos[random.pos.x,random.pos.y]+1
#svd
snp.svd.alter.2.add1.randpos <- svd(snp.genotype.alter.2.add1.randpos)
#comparison
(snp.svd.alter.2.add1.randpos$d-snp.svd$d)/snp.svd$d
#plot
pdf('snp.svd.alter.2.add1.randpos.pdf')
plot(((snp.svd.alter.2.add1.randpos$d-snp.svd$d)/snp.svd$d)[1:num_layers])
dev.off()

save(list=ls(),file='Change_Eigenvale.Rdata')

rand.add1.values <- array(0,c(10,127))
for (i in 1:10) {
# get a random position
random.pos.x <- random.pos(dim(snp.genotype)[1],length(main.pos.x))
random.pos.y <- random.pos(dim(snp.genotype)[2],length(main.pos.y))
# ## add 1
snp.genotype.alter.2.add1.randpos <- snp.genotype
snp.genotype.alter.2.add1.randpos[random.pos.x,random.pos.y] <- snp.genotype.alter.2.add1.randpos[random.pos.x,random.pos.y]+1
#svd
snp.svd.alter.2.add1.randpos <- svd(snp.genotype.alter.2.add1.randpos)
rand.add1.values[i,]<-snp.svd.alter.2.add1.randpos$d
print(i)
}

num_layers <- 10
probab <- array(0,c(num_layers,3))
for (layer in seq(num_layers)){
u <- filter.loadings(snp.svd$u[,layer],method='2',norm.cutoff=u.cutoff)
v <- filter.loadings(snp.svd$v[,layer],method='2',norm.cutoff=v.cutoff)
main.pos.x <- nonzero.pos(u)
main.pos.y <- nonzero.pos(v)
main.matrix <- snp.genotype[main.pos.x,main.pos.y]

count <- array(0,3)
names(count) <- c(0,1,2)
for (i in seq(dim(main.matrix)[1])){
    for (j in seq(dim(main.matrix)[2])){
        if (main.matrix[i,j]==0)
            count[1] <- count[1]+1
        else if (main.matrix[i,j]==1)
            count[2] <- count[2]+1
        else if (main.matrix[i,j]==2)
            count[3] <- count[3]+1
        else print(main.matrix[i,j])
    }
}
probab[layer,] <- count/sum(count)
print(layer)
print(probab[i])
}

# add 1 at the chosen positions
num_layers <- 10
change.sv <- array(0,c(num_layers,127))
for (layer in seq(num_layers)){
u <- filter.loadings(snp.svd$u[,layer],method='2',norm.cutoff=u.cutoff)
v <- filter.loadings(snp.svd$v[,layer],method='2',norm.cutoff=v.cutoff)
main.pos.x <- nonzero.pos(u)
main.pos.y <- nonzero.pos(v)
main.matrix <- snp.genotype[main.pos.x,main.pos.y]
change.matrix <- snp.genotype
change.matrix[main.pos.x,main.pos.y] <- change.matrix[main.pos.x,main.pos.y] + 1
change.sv[layer,] <- svd(change.matrix)$d
print(layer)
print(change.sv[layer,])
}

# add 1 at random positions
num_layers <- 10
rand.change.sv <- array(0,c(num_layers,127))
for (layer in seq(num_layers)){
u <- filter.loadings(snp.svd$u[,layer],method='2',norm.cutoff=u.cutoff)
v <- filter.loadings(snp.svd$v[,layer],method='2',norm.cutoff=v.cutoff)
main.pos.x <- nonzero.pos(u)
main.pos.y <- nonzero.pos(v)
random.pos.x <- random.pos(dim(snp.genotype)[1],length(main.pos.x))
random.pos.y <- random.pos(dim(snp.genotype)[2],length(main.pos.y))
rand.change.matrix <- snp.genotype
rand.change.matrix[random.pos.x,random.pos.y] <- rand.change.matrix[random.pos.x,random.pos.y] + 1
rand.change.sv[layer,] <- svd(rand.change.matrix)$d
print(layer)
print(rand.change.sv[layer,])
}

num.x.y <- length.pos.x * length.pos.y
num.scale <- mean(num.x.y[2:6])
num.x.y <- num.x.y/num.scale

change.rate <- change.sv
for (layer in seq(num_layers)){
    change.rate[layer,] <- ((change.sv[layer,]-snp.svd$d)/snp.svd$d)/num.x.y[layer]
}

rand.change.rate <- rand.change.sv
for (layer in seq(num_layers)){
    rand.change.rate[layer,] <- ((rand.change.sv[layer,]-snp.svd$d)/snp.svd$d)/num.x.y[layer]
}



t(change.sv[,1:10])
t(rand.change.sv[,1:10])
t(change.rate[,1:10])
t(rand.change.rate[,1:10])

library(pheatmap)
pdf('heatmap-change.pdf')
pheatmap(t(change.rate[1:6,1:6]),cluster_row=F,cluster_cols=F,border=F)
dev.off()
pdf('heatmap-rand.pdf')
pheatmap(t(rand.change.rate[1:6,1:6]),cluster_row=F,cluster_cols=F,border=F)
dev.off()

#####################
# subtract 1 
# subtract 1 at the chosen positions
num_layers <- 10
change.1.sv <- array(0,c(num_layers,127))
for (layer in seq(num_layers)){
u <- filter.loadings(snp.svd$u[,layer],method='2',norm.cutoff=u.cutoff)
v <- filter.loadings(snp.svd$v[,layer],method='2',norm.cutoff=v.cutoff)
main.pos.x <- nonzero.pos(u)
main.pos.y <- nonzero.pos(v)
main.matrix <- snp.genotype[main.pos.x,main.pos.y]
change.1.matrix <- snp.genotype
change.1.matrix[main.pos.x,main.pos.y] <- change.1.matrix[main.pos.x,main.pos.y] - 1
change.1.sv[layer,] <- svd(change.1.matrix)$d
print(layer)
print(change.1.sv[layer,])
}

# subtract 1 at random positions
num_layers <- 10
rand.change.1.sv <- array(0,c(num_layers,127))
for (layer in seq(num_layers)){
u <- filter.loadings(snp.svd$u[,layer],method='2',norm.cutoff=u.cutoff)
v <- filter.loadings(snp.svd$v[,layer],method='2',norm.cutoff=v.cutoff)
main.pos.x <- nonzero.pos(u)
main.pos.y <- nonzero.pos(v)
random.pos.x <- random.pos(dim(snp.genotype)[1],length(main.pos.x))
random.pos.y <- random.pos(dim(snp.genotype)[2],length(main.pos.y))
rand.change.1.matrix <- snp.genotype
rand.change.1.matrix[random.pos.x,random.pos.y] <- rand.change.1.matrix[random.pos.x,random.pos.y] - 1
rand.change.1.sv[layer,] <- svd(rand.change.1.matrix)$d
print(layer)
print(rand.change.1.sv[layer,])
}

change.1.rate <- change.1.sv
for (layer in seq(num_layers)){
    change.1.rate[layer,] <- ((change.1.sv[layer,]-snp.svd$d)/snp.svd$d)/num.x.y[layer]
}

rand.change.1.rate <- rand.change.1.sv
for (layer in seq(num_layers)){
    rand.change.1.rate[layer,] <- ((rand.change.1.sv[layer,]-snp.svd$d)/snp.svd$d)/num.x.y[layer]
}



t(change.1.sv[,1:10])
t(rand.change.1.sv[,1:10])
t(change.1.rate[,1:10])
t(rand.change.1.rate[,1:10])



library(pheatmap)
pdf('heatmap-change-1.pdf')
pheatmap(t(change.1.rate[1:6,1:6]),cluster_row=F,cluster_cols=F,border=F)
dev.off()
pdf('heatmap-rand-1.pdf')
pheatmap(t(rand.change.1.rate[1:6,1:6]),cluster_row=F,cluster_cols=F,border=F)
dev.off()
