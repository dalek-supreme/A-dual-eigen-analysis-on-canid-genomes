load('snp.svd.Rdata')
load('snp.genotype.Rdata')

# uniform noise
noise.level <- 0.005
noise.probability <- c(noise.level/2,noise.level/2,1-noise.level/2)
#generate uniform noise
rand.value <- array(runif(dim(snp.genotype)[1]*dim(snp.genotype)[2]),dim(snp.genotype))
noise.uniform <- array(-2,dim(snp.genotype))
for (i in seq(dim(snp.genotype)[1])){
    for (j in seq(dim(snp.genotype)[2])){
        if (rand.value[i,j] < noise.probability[1]) noise.uniform[i,j] <- 1
        else if (rand.value[i,j] < noise.probability[1]+noise.probability[2]) noise.uniform[i,j] <- -1
        else noise.uniform[i,j] <- 0
    }
}
#check
rand.count <- array(0,3)
for (i in seq(dim(noise.uniform)[1])){
    for (j in seq(dim(noise.uniform)[2])){
        if (noise.uniform[i,j]==1)
            rand.count[1] <- rand.count[1]+1
        else if (noise.uniform[i,j]==-1) rand.count[2] <- rand.count[2]+1
        else rand.count[3] <- rand.count[3]+1
    }
}
rand.probab <- rand.count/sum(rand.count)
(rand.probab - noise.probability)/noise.probability
#add noise
withnoise.genotype.uniform <- snp.genotype + noise.uniform
# svd
withnoise.svd.uniform <- svd(withnoise.genotype.uniform)

norm(noise.uniform,type='2')
# gaussian noise: seems unfeasible...
