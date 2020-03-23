library(gdsfmt)
library(SNPRelate)

vcf1.fn <- "Filtred_Published.vcf.bz2"
snpgdsVCF2GDS(vcf1.fn, "127_cnumref.gds", method = "copy.num.of.ref")

cat("Convertion Successful. Output: \"127_cnumref.gds\"")

vcf_ancient.fn <- "92indsAEDPCDTaimyr.biallelic.up.sort.vcf.gz"
snpgdsVCF2GDS(vcf_ancient.fn, "103_cnumref.gds", method = "copy.num.of.ref")

cat("Convertion Successful. Output: \"103_cnumref.gds\"")