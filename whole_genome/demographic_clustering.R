demographic.data <- read.csv('demographic.csv')
demographic.data.valid <- demographic.data[demographic.data$SampleID != 'Lcu2_Pastora',]

unique(demographic.data.valid$Breed)
unique(demographic.data.valid$Location)

S=cbind(demographic.data.valid$Longitude,demographic.data.valid$Latitude)
out.dist <- dist(S)
out.hclust <- hclust(out.dist)

num_clusters <- 5
out.id <- cutree(out.hclust,k=num_clusters)

cluster_geom <- list()
for(type in 1:num_clusters) cluster_geom[[type]]<-demographic.data.valid[out.id==type,]$SampleID

breeds <- unique(demographic.data.valid$Breed)
cluster_breed <- list()
for (type in seq_along(unique(demographic.data.valid$Breed)))
    cluster_breed [[type]] <- demographic.data.valid[demographic.data.valid$Breed==breeds[type],]$SampleID


cluster_dualeigen2 <- list()
k<-1
for (l1 in 2:10){
    for (l2 in 2:l1){
        for (l3 in 2:l2){
            for (l4 in 2:l3)
                cluster_dualeigen2[[k]] <- demographic.data.valid$SampleID[union(union(union(main.pos.x[[l1]],main.pos.x[[l2]]),main.pos.x[[l3]]),main.pos.x[[l4]])]
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


cluster_breed2 <- list()
k<-1
for (l1 in seq_along(breeds)){
    for (l2 in seq(l1)){
        for (l3 in seq(l2)){
            for (l4 in seq(l3))
                cluster_breed2[[k]] <- demographic.data.valid[demographic.data.valid$SampleID %in% union(union(union(cluster_breed[[l1]],cluster_breed[[l2]]),cluster_breed[[l3]]),cluster_breed[[l4]]),]$SampleID
                k <- k+1
        }
    }
}

cluster_index_breed <- list()
k<-1
for (l1 in seq_along(breeds)){
    for (l2 in seq(l1)){
        for (l3 in seq(l2)){
            for (l4 in seq(l3))
                cluster_index_breed[[k]] <- c(l1,l2,l3,l4)
                k <- k+1
        }
    }
}

##### cluster breed and cluster dualeigen comparison
identical_rate <- array(-1,c(length(cluster_breed2),length(cluster_dualeigen2)))
union_number <- array(-1,c(length(cluster_breed2),length(cluster_dualeigen2)))
intersect_number <- array(-1,c(length(cluster_breed2),length(cluster_dualeigen2)))
for (i in seq_along(cluster_breed2)){
    for (j in seq_along(cluster_dualeigen2)) {
        identical_rate[i,j] <- length(intersect(cluster_breed2[[i]],cluster_dualeigen2[[j]]))/length(union(cluster_breed2[[i]],cluster_dualeigen2[[j]]))
        union_number[i,j] <- length(union(cluster_breed2[[i]],cluster_dualeigen2[[j]]))
        intersect_number[i,j] <- length(intersect(cluster_breed2[[i]],cluster_dualeigen2[[j]]))
    }
}
max(identical_rate)

##### cluster breed and cluster tree comparison
identical_rate <- array(-1,c(length(cluster_breed2),length(cluster_tree)))
union_number <- array(-1,c(length(cluster_breed2),length(cluster_tree)))
intersect_number <- array(-1,c(length(cluster_breed2),length(cluster_tree)))
for (i in seq_along(cluster_breed2)){
    for (j in seq_along(cluster_tree)) {
        identical_rate[i,j] <- length(intersect(cluster_breed2[[i]],cluster_tree[[j]]))/length(union(cluster_breed2[[i]],cluster_tree[[j]]))
        union_number[i,j] <- length(union(cluster_breed2[[i]],cluster_tree[[j]]))
        intersect_number[i,j] <- length(intersect(cluster_breed2[[i]],cluster_tree[[j]]))
    }
}
max(identical_rate)


index_matrix <- array(list(),c(length(cluster_breed2),length(cluster_dualeigen2)))
for (i in seq_along(cluster_breed2)){
    for (j in seq_along(cluster_dualeigen2))
        index_matrix[i,j] <- list(c(i,j))
}

value_large <- identical_rate[identical_rate>0.9]

index_large <- index_matrix[identical_rate>0.9]

union_large <- union_number[identical_rate>0.9]

index_layer <- list()
for (i in seq_along(index_large)){
    index_layer[[i]] <- cluster_index[[index_large[[i]][2]]]
}

index_cluster <- list()
for (i in seq_along(index_large)){
    index_cluster[[i]] <- cluster_index_breed[[index_large[[i]][1]]]
}




breeds[index_cluster[[1]]]
main.pos.x[[3]]

length(demographic.data.valid[demographic.data.valid$Breed==breeds[index_cluster[[1]]][1],]$Breed)

demographic.data.valid[demographic.data.valid$Breed==breeds[index_cluster[[1]]][1],]$Breed