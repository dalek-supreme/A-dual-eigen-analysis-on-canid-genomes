prj <- function(vec_from,vec_to) {
    # obtains projection of vec_from on vec_to
    inner_prod <- sum(vec_from*vec_to)
    norm_vt <- sqrt(sum(vec_to*vec_to))
    return(inner_prod/norm_vt)
}


# count the SNP sites per gene
get.gene.snp.counts_and_sums <- function(layer,r_mtx_from_svd,gene.data) {
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
                #print(paste("in gene:",gene,"rank =",rank(r_mtx_from_svd[,layer])[site_num],sep=" "))
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

# take intersection with the genes with variations,
# as required by the wilcoxon enrichment analysis
get.gene.data.remove0s <- function(gene.stats,gene.data){
    zeros.index <-
        gene.data$ensembl_gene_id %in% 
            gene.stats[gene.stats$count==0,]$ensembl_gene_id
    return(gene.data[!zeros.index,])
}

#gene.stats[(gene.stats$ensembl_gene_id %in% unique(get.gene.data.remove0s(gene.stats,gene.data)$ensembl_gene_id)),]$count
#gene.stats[!(gene.stats$ensembl_gene_id %in% unique(get.gene.data.remove0s(gene.stats,gene.data)$ensembl_gene_id)),]$count

# removes zeros. takes average of loadings. and sorts.
get.gene.stats.sorted <- function(gene.stats){
    result.data <- data.frame(
        ensembl_gene_id=array("",length(gene.stats[gene.stats$count!=0,]$ensembl_gene_id)),
        avg_loading=array(0,length(gene.stats[gene.stats$count!=0,]$ensembl_gene_id)),
        stringsAsFactors = F
    )
    avg <- gene.stats[gene.stats$count!=0,]$sum/gene.stats[gene.stats$count!=0,]$count
    result.data$avg_loading <- sort(avg)
    result.data$ensembl_gene_id <- 
        gene.stats[gene.stats$count!=0,]$ensembl_gene_id[order(avg)]
    return(result.data)
}

# get subset for analysis with index from index.start to index.end
get.gene.subset <- function(gene.stats.sorted, index.start, index.end){
    return(gene.stats.sorted[index.start:index.end,])
}

#test:
#gene.subset<-get.gene.subset(get.gene.stats.sorted(gene.stats),1,100)

# does wilcoxon rank sum test on given subset w/ loadings
enrichment.wilcoxon <- function(gene.subset,gene.data){
    enrichment.result <- data.frame(
        go_id = array("",length(unique(gene.data$go_id))),
        p_value = array("",length(unique(gene.data$go_id))),
        stringsAsFactors = F
    )
    pathway.ids <- unique(gene.data$go_id)
    pathway.pval <- sapply(pathway.ids, function(x){
        index <- gene.data$go_id==x
        pathway.genes <- gene.data[index,]$ensembl_gene_id
        index <- gene.subset$ensembl_gene_id %in% pathway.genes
        loadings.in_pathway <- gene.subset[index,]$avg_loading
        loadings.out_of_pathway <- gene.subset[!index,]$avg_loading
        if (length(pathway.genes) >= 6){ # Question: why is this 6 and is this necessary and does it vary by different tests?
            wilcox.test(loadings.in_pathway, loadings.out_of_pathway, alternative = "greater")$p.val
        }
    })
}

enrichment.chi_square <- function()
enrichment.blahblah <- function()

# a generic function for enrichment
# does <enrichment.test> on given subset w/ loadings
# where <enrichment.test> is any statistical two-sample test,
# such as <wilcox.test> and the chi-square test, etc.
enrichment <- function(gene.subset,gene.data,enrichment.test){

}

# sort the enrichment result by p-value,
# then output the terms with the <num_terms> smallest p-values
enrichment.output <- function(enrichment.result, gene.data, num_terms){
    output <- data.frame(
        # information included in output
    )
}
