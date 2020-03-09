# enrichment analysis on loadings

setwd("~/enrichment") # set to the correct working dir

source("database.R")
source("utils.R")
#load("gene.db.Rdata")

sink("WILCOXON_enrichment.log")

# Wilcoxon Rank Sum Enrichment
layers <- 1:10
num_terms <- 50

for (layer in layers)
{
    gene.data <- gene.db.remove_unannotated(gene.db)

    gene.stats <- get.gene.snp.counts_and_sums(layer,s$v,gene.data)
    gene.data <- get.gene.data.remove0s(gene.stats,gene.data)
    
    gene.stats.sorted <-  get.gene.stats.sorted(gene.stats)
    gene.stats.positive <- gene.subset.pole_split(gene.stats.sorted,pole="positive")
    gene.stats.negative <- gene.subset.pole_split(gene.stats.sorted,pole="negative")

    result.positive <- enrichment.wilcoxon(gene.stats.positive,gene.data)
    result.negative <- enrichment.wilcoxon(gene.stats.negative,gene.data)

    setwd("~/enrichment/output/positive")
    filename <- paste("positive.wilcox.chr=",chr_t,".layer=",layer,".",sep="")
    enrichment.output(result.positive, gene.data, num_terms, filename, p_cutoff=0.05)
    setwd("~/enrichment/output/negative")   
    filename <- paste("negative.wilcox.chr=",chr_t,".layer=",layer,".",sep="")
    enrichment.output(result.negative, gene.data, num_terms, filename, p_cutoff=0.05)

    #enrichment.output.all(result, gene.data, filename)
    cat(paste("Layer",layer,"finished.\n",sep=" "))
}
sink()


