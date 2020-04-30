library(gdsfmt)
library(SNPRelate)

# for the 127-sample dataset first 
genome <- snpgdsOpen("127_biallelic.gds")

# # find which chromosome
# chromosomes <- read.gdsn(index.gdsn(genome, "snp.chromosome"))
# # find which sample
# samples <- read.gdsn(index.gdsn(genome, "sample.id"))
# # get the CanFam 3.1 coordinates
# snp.position <- read.gdsn(index.gdsn(genome, "snp.position"))

# missing_rate <- snpgdsSampMissRate(genome)


# snp.position

# # get indices at which chromosomes change
# ntypes<-1
# chr_type<-chromosomes[1]
# sizes<-1
# for (i in 2:length(snp.position)){
# 	if (chromosomes[i]!=chr_type[ntypes]){
# 		ntypes<-ntypes+1
# 		chr_type[ntypes]<-chromosomes[i]
# 		sizes[ntypes]<-i
# 	}
# }

# ntypes<-ntypes+1
# sizes[ntypes]<-i+1

# chromosomes[sizes]

# #SNP.index <- index.gdsn(genome, "genotype")
# SNP.matrix <- read.gdsn(index.gdsn(genome, "genotype"))

# # need to remove "3"s
# #???

# # LD pruning
# set.seed(1000)
# snpset <- snpgdsLDpruning(genome)

# SNP.svd <- svd(SNP.matrix)

# remove missing
snp.id.no_missing <- snpgdsSelectSNP(genome,missing.rate=0)
#SNP.matrix <- read.gdsn(index.gdsn(genome, "genotype"))
snp.genotype <- snpgdsGetGeno(genome,snp.id=snp.id.no_missing)
snp.genotype.pruned <- snpgdsGetGeno(genome,snp.id=snp.id.pruned)

snp.list.all <- snpgdsSNPList(genome)


#snp.list <- snp.list.all[snp.list.all$snp.id==snp.id.no_missing]

snp.list <- snp.list.all[snp.list.all$snp.id %in% snp.id.no_missing,]

snp.svd.test <- svd(snp.genotype[,1:1000])

snp.svd <- svd(snp.genotype)


#prune
load("snpset.Rdata")
snp.id.pruned <- unlist(unname(snpset))
snp.genotype.pruned <- snpgdsGetGeno(genome,snp.id=snp.id.pruned)
snp.list.pruned <- snp.list.all[snp.list.all$snp.id %in% snp.id.pruned,]

snp.svd.pruned <- svd(snp.genotype.pruned)


# to install
# source("http://bioconductor.org/biocLite.R")
# biocLite("pcaMethods")
# library(pcaMethods)
# library(matrixStats)

library(RobRSVD)

library(rsvd)

snp.svd.test.robust<-robustSvd(snp.genotype[,1:1000])

snp.svd.test.robust.reg <- RobRSVD(snp.genotype[,1:1000])

snp.svd.test.randsvd<-rsvd(snp.genotype[,1:1000])

snp.svd.randsvd <- rsvd(snp.genotype)
# Error in qr.default(Y, complete = FALSE) : too large a matrix for LINPACK

# # or
# snpid <- read.gdsn(index.gdsn(genome, "snp.id"))
# nmissid <- snpid[snpgdsSNPRateFreq(genome)$MissingRate == 0]

# grm <- snpgdsGRM(genome, snp.id=snpid, method="GCTA")
