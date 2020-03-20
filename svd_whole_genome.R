library(gdsfmt)
library(SNPRelate)

vcf1.fn <- "Filtred_Published.vcf.bz2"
snpgdsVCF2GDS(vcf1.fn, "127_cnumref.gds", method = "copy.num.of.ref")

cat("Convertion Successful. Output: \"127_cnumref.gds\"")

vcf_ancient.fn <- "92indsAEDPCDTaimyr.biallelic.up.sort.vcf.gz"
snpgdsVCF2GDS(vcf_ancient.fn, "103_cnumref.gds", method = "copy.num.of.ref")

cat("Convertion Successful. Output: \"103_cnumref.gds\"")

# for the 127-sample dataset first 
genome <- snpgdsOpen("127_biallelic.gds")

# find which chromosome
chromosomes <- read.gdsn(index.gdsn(genome, "snp.chromosome"))
# find which sample
samples <- read.gdsn(index.gdsn(genome, "sample.id"))
# get the CanFam 3.1 coordinates
snp.position <- read.gdsn(index.gdsn(genome, "snp.position"))

snp.position

# get indices at which chromosomes change
ntypes<-1
chr_type<-chromosomes[1]
sizes<-1
for (i in 2:length(snp.position)){
	if (chromosomes[i]!=chr_type[ntypes]){
		ntypes<-ntypes+1
		chr_type[ntypes]<-chromosomes[i]
		sizes[ntypes]<-i
	}
}

ntypes<-ntypes+1
sizes[ntypes]<-i+1

chromosomes[sizes]

#SNP.index <- index.gdsn(genome, "genotype")
SNP.matrix <- read.gdsn(index.gdsn(genome, "genotype"))

# need to remove "3"s
#???

SNP.svd <- svd(SNP.matrix)

