load('snp.genotype.Rdata')
# switch coding between 0 and 2
# compute svd using sparsesvd() 
# ...

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

save(snp.genotype.recoding,file='snp.genotype.recoding.Rdata')

snp.recoding.svd <- svd(snp.genotype.recoding)
save(snp.recoding.svd,file='snp.recoding.svd.Rdata')
