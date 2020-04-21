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
# layer 2
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




