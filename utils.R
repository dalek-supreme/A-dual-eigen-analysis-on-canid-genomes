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
                if(site_num>dim(r_mtx_from_svd)[1]){
                    break
                }
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
                if(site_num>dim(r_mtx_from_svd)[1]){
                    break
                }
            }
        if(site_num>dim(r_mtx_from_svd)[1]){
            break
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


<<<<<<< HEAD
<<<<<<< HEAD
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
        print(paste("GO term:",x,sep=""))
        if (length(pathway.genes) >= 6 # Question: why is this 6 and is this necessary and does it vary by different tests?
            && length(loadings.in_pathway)>0
            && length(loadings.out_of_pathway)>0
        ){ 
            if(NA %in% c(loadings.in_pathway,loadings.out_of_pathway)){
                print("NA!!!!")
            }
            print(paste("GO term:",x," beginning wilcox.test",sep=""))
            print(index)
            print(loadings.in_pathway)
            print(loadings.out_of_pathway)
            print(wilcox.test(
                as.numeric(loadings.in_pathway), 
                as.numeric(loadings.out_of_pathway), 
                alternative = "greater"
            )$p.val)
            wilcox.test(
                as.numeric(loadings.in_pathway), 
                as.numeric(loadings.out_of_pathway), 
                alternative = "greater"
                )$p.val
        }
    })
    enrichment.result$go_id <- pathway.ids
    enrichment.result$p_value <- pathway.pval
    return(enrichment.result)
}
=======
>>>>>>> parent of c5e34e9... Update utils.R
=======
>>>>>>> parent of c5e34e9... Update utils.R

pathway.id <- unique(GO.data[,2])


enrichment.wilcoxon <- function(experiment.data, pathway.data) {

<<<<<<< HEAD
<<<<<<< HEAD
# sort the enrichment result by p-value,
# then output the terms with the <num_terms> smallest p-values
enrichment.output <- function(enrichment.result, gene.data, num_terms, filename){
    output.data <- data.frame(
        # information included in output
        p_value=array("",num_terms),
        go_id=array("",num_terms),
        name_1006=array("",num_terms),
        go_linkage_type=array("",num_terms),
        namespace_1003=array("",num_terms),
        definition_1006=array("",num_terms),
        stringsAsFactors = F
    )
    pval <- enrichment.result$p_value
    output.data$p_value <- sort(pval)[1:num_terms]
    output.data$go_id <- enrichment.result$go_id[order(pval)[1:num_terms]]
    for (id in output.data$go_id){
        output.data$name_1006 <- unique(gene.data[gene.data$go_id==id,]$name_1006)
        output.data$go_linkage_type <- unique(gene.data[gene.data$go_id==id,]$go_linkage_type)
        output.data$namespace_1003 <- unique(gene.data[gene.data$go_id==id,]$namespace_1003)
        output.data$definition_1006 <- unique(gene.data[gene.data$go_id==id,]$definition_1006)
    }
    write.csv(output.data,file=filename)
    return(output.data)
}
=======
=======
>>>>>>> parent of c5e34e9... Update utils.R
    pathway.pval <- sapply(pathway.id, function(x){
        
    })
    # result<-data.frame(...)
    # return(result)
<<<<<<< HEAD
}
>>>>>>> parent of c5e34e9... Update utils.R
=======
}
>>>>>>> parent of c5e34e9... Update utils.R
