cluster_from_tree <- read.csv('cluster_results',sep='\t')
demographic.data <- read.csv('demographic.csv')
# cluster_from_tree[cluster_from_tree$ClusterNumber==5,]$SequenceName

## Get main positions:
u.cutoff <- 0.75 -> v.cutoff
main.pos.x <- list()
main.pos.y <- list()
for (layer in seq(num_layers)) {
    u <- filter.loadings(snp.svd$u[,layer],method='2',norm.cutoff=u.cutoff)
    v <- filter.loadings(snp.svd$v[,layer],method='2',norm.cutoff=v.cutoff)
    main.pos.x[[layer]] <- nonzero.pos(u)
    main.pos.y[[layer]] <- nonzero.pos(v)
}

cluster_tree<-list()
for(type in 1:5) cluster_tree[[type]]<-cluster_from_tree[cluster_from_tree$ClusterNumber==type,]$SequenceName

cluster_dualeigen<-list()
for (layer in 2:10) cluster_dualeigen[[layer]] <- demographic.data$SampleID[main.pos.x[[layer]]]

identical_rate <- array(-1,c(length(cluster_tree),length(cluster_dualeigen)))
for (i in seq_along(cluster_tree)){
    for (j in seq_along(cluster_dualeigen))
        identical_rate[i,j] <- length(intersect(cluster_tree[[i]],cluster_dualeigen[[j]]))/length(union(cluster_tree[[i]],cluster_dualeigen[[j]]))
}

save(identical_rate_0.8,identical_rate_0.9,identical_rate_0.85,identical_rate_0.7,identical_rate_0.75,file='cluster_comparison.Rdata')
