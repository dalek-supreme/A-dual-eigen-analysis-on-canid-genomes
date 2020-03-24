# enrichment analysis on loadings

source("database.R")
source("utils.R")
#load("gene.db.Rdata")

#sink("WILCOXON_enrichment.log")

# Wilcoxon Rank Sum Enrichment
num_layers <- 10
layers <- seq(num_layers)
num_terms <- 50

system(paste("NUM_LAYERS=",num_layers,sep=""))
system("bash create_result_dir.sh")

t1 <- proc.time()

gene.stats <- gene.loading.stats(snp.svd$v,gene.db)

# load("gene.db.Rdata")
# load("gene.stats.Rdata")
t5<-proc.time()
for (layer in layers)
{
    t3 <- proc.time()
    gene.stats.layer <- gene.loadings.stats.get_layer(gene.stats,layer)
    gene.stats.sorted <-  gene.stats.sort(gene.stats)
    gene.stats.positive <- gene.subset.pole_split(gene.stats.sorted,pole="positive")
    gene.stats.negative <- gene.subset.pole_split(gene.stats.sorted,pole="negative")

    result.positive <- enrichment.wilcoxon(gene.stats.positive,gene.db) # about 7 minutes
    result.negative <- enrichment.wilcoxon(gene.stats.negative,gene.db)
    setwd(paste("~/enrichment/",layer,sep=""))
    save(result.positive,result.negative,file="result.Rdata")
    setwd(paste("~/enrichment/",layer,"positive",sep=""))
    enrichment.output(result.positive, gene.db, num_terms, filename="wilcox", p_cutoff=0.05)
    setwd(paste("~/enrichment/",layer,"negative",sep=""))
    enrichment.output(result.negative, gene.db, num_terms, filename="wilcox", p_cutoff=0.05)
    #enrichment.output.all(result, gene.data, filename)
    t4 <- proc.time()-t3
    print(paste("Layer ",layer," finished. Time elapsed: ",t4[3],sep=""))
    t6<-proc.time()-t5
    print(paste("Estimated time remaining: ",t6[3]*(length(layers)/layer-1)/60,"minutes.",sep=""))
}

t7 <- proc.time()-t0
print(paste("Analysis complete. Total time: ",t7[3],sep=""))

