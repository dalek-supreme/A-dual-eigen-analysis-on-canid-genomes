# # LD pruning
set.seed(1000)
snpset <- snpgdsLDpruning(genome)
save(snpset, file="snpset.Rdata")
cat("LD pruning completed. data saved.")
snp.id.pruned <- unlist(unname(snpset))