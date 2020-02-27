# enrichment analysis on loadings
source(utils.R)
source(database.R)

#setwd("") # set to the correct working dir

load("gene.data.Rdata")

layers <- 1:10
num_terms <- 10
filename <- paste("output_wilcox_chr=",chr_t,".csv",sep="")
t1<-Sys.time()
for (layer in layers){
    gene.stats <- get.gene.snp.counts_and_sums(layer,s$v,gene.data)
    gene.data <- get.gene.data.remove0s(gene.stats,gene.data)
    gene.stats.sorted <-  get.gene.stats.sorted(gene.stats)
    gene.subset <- get.gene.subset(gene.stats.sorted,1,100)
    result <- enrichment.wilcoxon(gene.subset,gene.data)
    filename <- paste("output_wilcox_chr=",chr_t,"_layer=",layer,".csv",sep="")
    output <- enrichment.output(result, gene.data, num_terms, filename, p_cutoff=0.05)
}
t2<-Sys.time()
print(paste("Time elapsed: ",t2-t1,sep=""))