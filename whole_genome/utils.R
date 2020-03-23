gene.loading.stats <- function(snp.svd.right, gene.data) {
    gene.isunique <- !duplicated(gene.data$ensembl_gene_id)
    gene.unique.index <- seq_along(gene.data$ensembl_gene_id)[gene.isunique]
    gene.unique <- data.frame(
                    ensembl_gene_id = gene.data$ensembl_gene_id[gene.unique.index],
                    chromosome_name=gene.data$chromosome_name[gene.unique.index],
                    start_position = gene.data$start_position[gene.unique.index],
                    end_position = gene.data$end_position[gene.unique.index],
                    stringsAsFactors = F
                )
    
    result <- data.frame(
        ensembl_gene_id = gene.unique$ensembl_gene_id,
        loading.count = array(list(array(0,dim(snp.svd.right)[2])),length(gene.unique$ensembl_gene_id)),
        loading.sum = array(list(array(0,dim(snp.svd.right)[2])),length(gene.unique$ensembl_gene_id)),
        loading.avg = array(list(array(0,dim(snp.svd.right)[2])),length(gene.unique$ensembl_gene_id)),
        stringsAsFactors = F
    )

    for (gene.index in seq_along(gene.unique$ensembl_gene_id)) {
        gene.loading.count <- 0
        gene.loading.sum <- array(0,dim(snp.svd.right)[2]) # do it once for all layers
        gene.loading.avg <- array(0,dim(snp.svd.right)[2])
        for (site.index in dim(snp.svd.right)[1]) {
            if (snp.list$chromosome[site.index] == gene.unique$chromosome_name[gene.index]){
                if (gene.unique$start_position[gene.index] <= snp.list$position[site.index]
                && gene.unique$end_position[gene.index] >= snp.list$position[site.index]
                ) {
                    gene.loading.count <- gene.loading.count+1
                    gene.loading.sum <- gene.loading.sum + snp.svd.right[site.index,]
                }
            }
        }
        if (gene.loading.count != 0) {
            gene.loading.avg <- gene.loading.sum/gene.loading.count
        }
        result$loading.count[[gene.index]] <- list(gene.loading.count)
        result$loading.sum[[gene.index]] <- list(gene.loading.sum)
        result$loading.avg[[gene.index]] <- list(gene.loading.avg)
    }
}
