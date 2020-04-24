load('snp.svd.Rdata')
load('snp.genotype.Rdata')
repeat_time <- 10
noise.norm <- array(0,repeat_time)
noise.sv <- array(0,c(repeat_time,127))
for (num_repeat in seq(repeat_time)){
# uniform noise
noise.level <- 0.01
noise.probability <- c(noise.level/2,noise.level/2,1-noise.level)
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

noise.norm[num_repeat]<-norm(noise.uniform,type='2')
noise.sv[num_repeat,]<-withnoise.svd.uniform$d
# gaussian noise: seems unfeasible...
print(num_repeat)
print(noise.norm)
}
save(noise.level,noise.norm,noise.sv,file='noise11.Rdata')

# for (i in 1:10){
#     load(paste('noise',i,sep=''))
#     print(paste(noise.level/2,mean(noise.norm),mean(noise.sv[,5],)))
# }

unit.vec.angle <- function(v1,v2,unit=c('rad','deg')){
    if (unit == 'rad') return(acos(sum(v1*v2)))
    if (unit == 'deg') return(acos(sum(v1*v2))/pi*180)
    else stop('Wrong argument for unit: use \'rad\' or \'deg\'\n')
}

angle.between <- array(-1,c(10,2))
for (layer in 1:10){
    angle.between[layer,1] <- unit.vec.angle(snp.svd$u[,layer],withnoise.svd.uniform$u[,layer],'deg')
    angle.between[layer,2] <- unit.vec.angle(snp.svd$v[,layer],withnoise.svd.uniform$v[,layer],'deg')
    print(layer)
}

             [,1]       [,2]
 [1,] 0.004419025  0.1861792
 [2,] 0.222643350  7.1497291
 [3,] 0.265796808  8.0103817
 [4,] 0.402121451  9.7126283
 [5,] 0.678335774 11.0517146
 [6,] 0.687942768 11.1871505
 [7,] 0.862722169 12.3831489
 [8,] 1.202461808 12.8328358
 [9,] 1.060961290 13.2424848
[10,] 1.183355535 13.6761971

             [,1]       [,2]
 [1,] 0.005414321  0.2623535
 [2,] 0.296407771 10.0515599
 [3,] 0.379023795 11.1997354
 [4,] 0.625800599 13.6036375
 [5,] 1.275941254 15.3226323
 [6,] 1.215029737 15.5223180
 [7,] 1.062603895 17.2273549
 [8,] 2.097589886 17.8481633
 [9,] 2.455247767 18.5656865
[10,] 2.067383514 18.8647519