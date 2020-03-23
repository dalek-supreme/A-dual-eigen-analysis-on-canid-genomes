library(gdsfmt)
library(SNPRelate)

# for the 127-sample dataset first 
genome <- snpgdsOpen("127_biallelic.gds")
snp.list.all <- snpgdsSNPList(genome)

snp.id.no_missing <- snpgdsSelectSNP(genome,missing.rate=0)

load("snpset.Rdata")
snp.id.pruned <- unlist(unname(snpset))
snp.id.pruned.no_missing <- intersect(snp.id.pruned,snp.id.no_missing)

snp.genotype <- snpgdsGetGeno(genome,snp.id=snp.id.pruned.no_missing)
snp.list <- snp.list.all[snp.list.all$snp.id %in% snp.id.pruned.no_missing,]

snp.svd <- svd(snp.genotype)

