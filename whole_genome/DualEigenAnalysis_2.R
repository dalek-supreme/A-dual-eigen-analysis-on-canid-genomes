load('gene.db.Rdata')
source('utils.R')
source('DualEigen_utils.R')
load('snp.svd.Rdata')
load('snp.genotype.Rdata')

#sink("WILCOXON_enrichment.log")

# Wilcoxon Rank Sum Enrichment
num_layers <- 10
layers <- 2:10
num_terms <- 50

system(paste("NUM_LAYERS=",num_layers,sep=""))
system("bash create_result_dir.sh")

t1 <- proc.time()

gene.stats <- gene.loading.stats(snp.svd$v,gene.db)
########## stuck here ##########
save(gene.stats,file='gene.stats.Rdata')

## Get main positions:
u.cutoff <- 0.8 -> v.cutoff
main.pos.x <- list()
main.pos.y <- list()
for (layer in seq(num_layers)) {
    u <- filter.loadings(snp.svd$u[,layer],method='2',norm.cutoff=u.cutoff)
    v <- filter.loadings(snp.svd$v[,layer],method='2',norm.cutoff=v.cutoff)
    main.pos.x[[layer]] <- nonzero.pos(u)
    main.pos.y[[layer]] <- nonzero.pos(v)
}


# load("gene.db.Rdata")
# load("gene.stats.Rdata")
t5<-proc.time()
for (layer in layers)
{
    t3 <- proc.time()
    gene.stats.layer <- gene.loadings.stats.get_layer(gene.stats,layer)
    gene.stats.sorted <-  gene.stats.sort(gene.stats.layer)
    # gene.stats.positive <- gene.subset.pole_split(gene.stats.sorted,pole="positive")
    # gene.stats.negative <- gene.subset.pole_split(gene.stats.sorted,pole="negative")
    gene.stats.main <- gene.subset.select(gene.stats.sorted, main.pos.y[[layer]])
    result.main <- enrichment.wilcoxon(gene.stats.main,gene.db) # about 7 minutes
    setwd(paste("~/enrichment2/",layer,sep=""))
    save(result.main,,file="result.Rdata")
    setwd(paste("~/enrichment2/",layer,sep=""))
    enrichment.output(result.main, gene.db, num_terms, filename="wilcox", p_cutoff=0.05)
    t4 <- proc.time()-t3
    print(paste("Layer ",layer," finished. Time elapsed: ",t4[3],sep=""))
    t6<-proc.time()-t5
    print(paste("Estimated time remaining: ",t6[3]*(length(layers)/layer-1)/60,"minutes.",sep=""))
}

t7 <- proc.time()-t0
print(paste("Analysis complete. Total time: ",t7[3],sep=""))

