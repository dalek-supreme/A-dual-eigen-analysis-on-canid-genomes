load('snp.genotype.Rdata')
load('snp.svd.Rdata')
source('DualEigen_utils.R')
u.cutoff <- 0.8 -> v.cutoff
num_layers <- 8
B <- snp.svd$d[1]*snp.svd$u[,1]%*%t(snp.svd$v[,1])
for (layer in 2:num_layers){
    u <- filter.loadings(snp.svd$u[,layer],method='2',norm.cutoff=u.cutoff)
    v <- filter.loadings(snp.svd$v[,layer],method='2',norm.cutoff=v.cutoff)
    B <- B + snp.svd$d[layer]*u %*% t(v)
}
err <- norm(snp.genotype-B, type='2')



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

