prj <- function(vec_from,vec_to) {
    # obtains projection of vec_from on vec_to
    inner_prod <- sum(vec_from*vec_to)
    norm_vt <- sqrt(sum(vec_to*vec_to))
    return(inner_prod/norm_vt)
}


# count the SNP sites per gene
gene.stats <- function(layer,r_mtx_from_svd,gene.data) {
    result.data <- data.frame(
        ensembl_gene_id = unique(gene.data$ensembl_gene_id),
        count = array(0,length(unique(gene.data$ensembl_gene_id))),
        sum = array(0,length(unique(gene.data$ensembl_gene_id))),
        stringsAsFactors = F
    )
    site_num<-1
    for (gene in unique(gene.data$ensembl_gene_id)){
        snp.loading.count <- 0
        snp.loading.sum <- 0
        # locate site_num pointer into the gene, break if no
        while (
            pos[coordinates[site_num]] < 
            unique(gene.data[gene.data$ensembl_gene_id==gene,]$start_position)
            ){
                site_num <- site_num+1
                if (
                    pos[coordinates[site_num]] >
                    unique(gene.data[gene.data$ensembl_gene_id==gene,]$end_position)
                ){
                    #print(paste("No SNP sites in gene ",gene,sep=""))
                    break
                }
            }
        while (
            pos[coordinates[site_num]] <=
            unique(gene.data[gene.data$ensembl_gene_id==gene,]$end_position)
            ){
                #print(paste("found in gene: ",gene," at site: ",site_num,". ", site_num/dim(r_mtx_from_svd)[1]*100,"% completed.",sep=""))
                snp.loading.count <- snp.loading.count + 1
                snp.loading.sum <- snp.loading.sum + r_mtx_from_svd[site_num,layer]
                site_num <- site_num+1
            }
        result.data[result.data$ensembl_gene_id==gene,]$count <- snp.loading.count
        result.data[result.data$ensembl_gene_id==gene,]$sum <- snp.loading.sum
    }
    # ignore this error message when the loop is terminated as the iteration on site_num is complete:
    # "Error in while (pos[coordinates[site_num]] <= unique(gene.data[gene.data$ensembl_gene_id ==  : 
    #  missing value where TRUE/FALSE needed"
    return(result.data)
}

enrichment.wilcoxon <- function()
enrichment.chi_square <- function()
enrichment.blahblah <- function()



pathway.id <- unique(GO.data[,2])


enrichment.wilcoxon <- function(experiment.data, pathway.data) {

    pathway.pval <- sapply(pathway.id, function(x){
        
    })
    # result<-data.frame(...)
    # return(result)
}