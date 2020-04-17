library(gdsfmt)
library(SNPRelate)

# for the 127-sample dataset first 
genome <- snpgdsOpen("127_biallelic.gds")
sample.id <- read.gdsn(index.gdsn(genome, "sample.id"))

snp.list.all <- snpgdsSNPList(genome)

snp.id.no_missing <- snpgdsSelectSNP(genome,missing.rate=0)

load("snpset.Rdata")
#source("LD_pruning.R")
snp.id.pruned <- unlist(unname(snpset))
snp.id.pruned.no_missing <- intersect(snp.id.pruned,snp.id.no_missing)

snp.genotype <- snpgdsGetGeno(genome,snp.id=snp.id.pruned.no_missing)
snp.list <- snp.list.all[snp.list.all$snp.id %in% snp.id.pruned.no_missing,]

save(snp.genotype,file='snp.genotype.Rdata')

snp.svd <- svd(snp.genotype)
save(snp.svd,file='snp.svd.Rdata')


# randomized svd
library(rsvd)

snp.rsvd <- rsvd(snp.genotype)
save(snp.rsvd,file='snp.rsvd.Rdata')

# robust svd
library(pcaMethods)
snp.robsvd <- robustSvd(snp.genotype)
save(snp.robsvd,file='snp.robsvd.Rdata')
# memory exceeded!!!

library(RobRSVD)
snp.RobRSVD1 <- RobRSVD(snp.genotype[1:10000],irobust=TRUE, istablize=F)
save(snp.RobRSVD,file='snp.RobRSVD.Rdata')
# doesnt'work!!!