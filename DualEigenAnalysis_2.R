############################
library(gdsfmt)
library(SNPRelate)

# for the 127-sample dataset first 
genome <- snpgdsOpen("127_biallelic.gds")

snp.list.all <- snpgdsSNPList(genome)

snp.id.no_missing <- snpgdsSelectSNP(genome,missing.rate=0)

load("snpset.Rdata")
#source("LD_pruning.R")
snp.id.pruned <- unlist(unname(snpset))
snp.id.pruned.no_missing <- intersect(snp.id.pruned,snp.id.no_missing)

snp.genotype <- snpgdsGetGeno(genome,snp.id=snp.id.pruned.no_missing)
snp.list <- snp.list.all[snp.list.all$snp.id %in% snp.id.pruned.no_missing,]

snp.svd <- svd(snp.genotype)

sample.id <- read.gdsn(index.gdsn(genome, "sample.id"))

print('svd finished.')

############################

load('gene.db.Rdata')
source('utils-2.R')
source('DualEigen_utils.R')
#load('snp.svd.Rdata')
#load('snp.genotype.Rdata')

#sink("WILCOXON_enrichment.log")

# Wilcoxon Rank Sum Enrichment
num_layers <- 10
layers <- 2:10
num_terms <- 50

# system(paste("NUM_LAYERS=",num_layers,sep=""))
# system("bash create_result_dir.sh")

t0 <- proc.time()
t1 <- proc.time()

# gene.stats <- gene.loading.stats(snp.svd$v,gene.db)
# ########## stuck here ##########
# save(gene.stats,file='gene.stats2.Rdata')

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

for (layer in seq(num_layers)) {
    print(paste(layer,length(main.pos.x[[layer]]),length(main.pos.y[[layer]])))
}

print('main positions found.')

# load("gene.db.Rdata")

t5<-proc.time()
for (layer in layers)
{
    t3 <- proc.time()
    gene.stats <- gene.loading.stats.1_layer_only(snp.svd$v[main.pos.y[[layer]],layer],gene.db)
    
    #gene.stats.layer <- gene.loadings.stats.get_layer(gene.stats)
    #gene.stats.sorted <-  gene.stats.sort(gene.stats.layer)
    # gene.stats.positive <- gene.subset.pole_split(gene.stats.sorted,pole="positive")
    # gene.stats.negative <- gene.subset.pole_split(gene.stats.sorted,pole="negative")
    #gene.stats.main <- gene.subset.select(gene.stats.layer, main.pos.y[[layer]])
    result.molecular_function <- enrichment.wilcoxon(gene.stats,gene.db,'molecular_function') # about 7 minutes
    result.biological_process <- enrichment.wilcoxon(gene.stats,gene.db,'biological_process') # about 7 minutes
    result.cellular_component <- enrichment.wilcoxon(gene.stats,gene.db,'cellular_component') # about 7 minutes
    
    setwd("~/enrichment5/")
    save(result.molecular_function,result.biological_process,result.cellular_component,file=paste('result-',layer,'.Rdata',sep=''))
    enrichment.output(result.molecular_function, gene.db, num_terms, filename=paste('wilcox-molecular_function',layer,sep=''), p_cutoff=0.05)
    enrichment.output(result.biological_process, gene.db, num_terms, filename=paste('wilcox-biological_process',layer,sep=''), p_cutoff=0.05)
    enrichment.output(result.cellular_component, gene.db, num_terms, filename=paste('wilcox-cellular_component',layer,sep=''), p_cutoff=0.05)
    


    t4 <- proc.time()-t3
    print(paste("Layer ",layer," finished. Time elapsed: ",t4[3],sep=""))
    t6<-proc.time()-t5
    print(paste("Estimated time remaining: ",t6[3]*(length(layers)/layer-1)/60,"minutes.",sep=""))
}

t7 <- proc.time()-t0
print(paste("Analysis complete. Total time: ",t7[3],sep=""))

