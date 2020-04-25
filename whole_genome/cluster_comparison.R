cluster_from_tree <- read.csv('cluster_results',sep='\t')
demographic.data <- read.csv('demographic.csv')
# cluster_from_tree[cluster_from_tree$ClusterNumber==5,]$SequenceName

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

cluster_tree<-list()
for(type in 1:5) cluster_tree[[type]]<-cluster_from_tree[cluster_from_tree$ClusterNumber==type,]$SequenceName

cluster_dualeigen<-list()
for (layer in 2:10) cluster_dualeigen[[layer]] <- demographic.data$SampleID[main.pos.x[[layer]]]

cluster_dualeigen2 <- list()
k<-1
for (l1 in 2:10){
    for (l2 in 2:l1){
        for (l3 in 2:l2){
            for (l4 in 2:l3)
                cluster_dualeigen2[[k]] <- demographic.data$SampleID[union(union(union(main.pos.x[[l1]],main.pos.x[[l2]]),main.pos.x[[l3]]),main.pos.x[[l4]])]
                k <- k+1
        }
    }
}

cluster_index <- list()
k<-1
for (l1 in 2:10){
    for (l2 in 2:l1){
        for (l3 in 2:l2){
            for (l4 in 2:l3)
                cluster_index[[k]] <- c(l1,l2,l3,l4)
                k <- k+1
        }
    }
}

identical_rate <- array(-1,c(length(cluster_tree),length(cluster_dualeigen2)))
for (i in seq_along(cluster_tree)){
    for (j in seq_along(cluster_dualeigen2))
        identical_rate[i,j] <- length(intersect(cluster_tree[[i]],cluster_dualeigen2[[j]]))/length(union(cluster_tree[[i]],cluster_dualeigen2[[j]]))
}
max(identical_rate)

index_matrix <- array(list(),c(length(cluster_tree),length(cluster_dualeigen2)))
for (i in seq_along(cluster_tree)){
    for (j in seq_along(cluster_dualeigen2))
        index_matrix[i,j] <- list(c(i,j))
}

value_large <- identical_rate[identical_rate>0.7]

index_large <- index_matrix[identical_rate>0.7]

index_layer <- list()
for (i in seq_along(index_large)){
    index_layer[[i]] <- cluster_index[[index_large[[i]][2]]]
}

index_cluster <- array()
for (i in seq_along(index_large)){
    index_cluster[i] <- index_large[[i]][1]
}


save(identical_rate_0.8,identical_rate_0.9,identical_rate_0.85,identical_rate_0.7,identical_rate_0.75,file='cluster_comparison.Rdata')
