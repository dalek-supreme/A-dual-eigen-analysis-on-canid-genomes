count <- array(0,3)
names(count) <- c(0,1,2)
for (i in seq(dim(snp.genotype)[1])){
    for (j in seq(dim(snp.genotype)[2])){
        if (snp.genotype[i,j]==0)
            count[1] <- count[1]+1
        else if (snp.genotype[i,j]==1)
            count[2] <- count[2]+1
        else if (snp.genotype[i,j]==2)
            count[3] <- count[3]+1
        else print(snp.genotype[i,j])
    }
}
probab <- count/sum(count)
rand.value <- array(runif(dim(snp.genotype)[1]*dim(snp.genotype)[2]),dim(snp.genotype))
rand.genotype <- array(0,dim(snp.genotype))
for (i in dim(snp.genotype)[1]){
    for (j in dim(snp.genotype)[2]){
        if (rand.value[i][j] < probab[1]) rand.genotype[i][j] <- 0
        else if (rand.value[i][j] < prabab[2]) rand.genotype[i][j] <- 1
        else rand.genotype[i][j] <- 2
    }
}

rand.svd <- svd(rand.genotype)