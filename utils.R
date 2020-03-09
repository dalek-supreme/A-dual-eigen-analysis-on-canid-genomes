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
                #print(paste("in gene:",gene,"rank =",rank(r_mtx_from_svd[,layer])[site_num],sep=" "))
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

# take intersection with the genes with variations,
# as required by the wilcoxon enrichment analysis
get.gene.data.remove0s <- function(gene.stats,gene.data){
    zeros.index <-
        gene.data$ensembl_gene_id %in% 
            gene.stats[gene.stats$count==0,]$ensembl_gene_id
    return(gene.data[!zeros.index,])
}

# remove genes with no GO ids
gene.db.remove_unannotated <- function(gene.data){
    return(gene.data[gene.data$go_id!="",])
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
get.gene.subset <- function(gene.stats.sorted, index.range){
    if (NA %in% gene.stats.sorted[index.range,1]){
        stop(paste("get.gene.subset: index.start or index.end out of bound! Legit range: ",range(gene.stats.sorted[,1]),sep=""))
    }
    return(gene.stats.sorted[index.range,])
}

gene.subset.pole_split <- function(gene.stats,pole){
    if (pole=="positive"){
        return(gene.stats[gene.stats$avg_loading>0,])
    } else if (pole=="negative") {
        return(gene.stats[gene.stats$avg_loading<0,])
    } else {
        stop(paste(" Invalid argument \"",pole,"\": should be \"positive\" or \"negative\"!\n",sep=""))
    }
}

#test:
#gene.subset<-get.gene.subset(get.gene.stats.sorted(gene.stats),1,100)

# does wilcoxon rank sum test on given subset w/ loadings
enrichment.wilcoxon <- function(gene.subset,gene.data){
    enrichment.result <- data.frame(
        go_id = array("",length(unique(gene.data$go_id))),
        p_value.greater = array("",length(unique(gene.data$go_id))),
        p_value.less = array("",length(unique(gene.data$go_id))),
        p_value.two.sided = array("",length(unique(gene.data$go_id))),
        genes.in.set.count = array("",length(unique(gene.data$go_id))),
        genes.out_of.set.count = array("",length(unique(gene.data$go_id))),
        stringsAsFactors = F
    )
    pathway.ids <- unique(gene.data$go_id)

    for (id_num in (1:length(pathway.ids))){
        index <- (gene.data$go_id==pathway.ids[id_num])
        pathway.genes <- gene.data[index,]$ensembl_gene_id
        index <- gene.subset$ensembl_gene_id %in% pathway.genes
        loadings.in_pathway <- gene.subset[index,]$avg_loading
        loadings.out_of_pathway <- gene.subset[!index,]$avg_loading
        genes.in.set.count <- length(loadings.in_pathway)
        genes.out_of.set.count <- length(loadings.out_of_pathway)

        enrichment.result$genes.in.set.count[id_num] <- genes.in.set.count
        enrichment.result$genes.out_of.set.count[id_num] <- genes.out_of.set.count

        if (length(pathway.genes) >= 6 # Question: why is this 6 and is this necessary and does it vary by different tests?
            && genes.in.set.count>0
            && genes.out_of.set.count>0
        ){ 
            enrichment.result$p_value.greater[id_num] <-    
                wilcox.test(
                as.numeric(loadings.in_pathway), 
                as.numeric(loadings.out_of_pathway), 
                alternative = "greater"
                )$p.val
            enrichment.result$p_value.less[id_num] <-    
                wilcox.test(
                as.numeric(loadings.in_pathway), 
                as.numeric(loadings.out_of_pathway), 
                alternative = "less"
                )$p.val
            enrichment.result$p_value.two.sided[id_num] <-    
                wilcox.test(
                as.numeric(loadings.in_pathway), 
                as.numeric(loadings.out_of_pathway), 
                alternative = "two.sided"
                )$p.val
        } else {
            enrichment.result$p_value.greater[id_num] <- 10 # set as 10>1 so as to be omitted after sorting 
            enrichment.result$p_value.less[id_num] <- 10
            enrichment.result$p_value.two.sided[id_num] <- 10
        }
    }

    enrichment.result$go_id <- pathway.ids
    return(enrichment.result)
}

# # a generic function for enrichment
# # does <enrichment.test> on given subset w/ loadings
# # where <enrichment.test> is any statistical two-sample test,
# # such as <wilcox.test> and the chi-square test, etc.
# enrichment <- function(gene.subset,gene.data,stat.test){
#     enrichment.result <- data.frame(
#         go_id = array("",length(unique(gene.data$go_id))),
#         p_value = array("",length(unique(gene.data$go_id))),
#         stringsAsFactors = F
#     )
#     pathway.ids <- unique(gene.data$go_id)
#     pathway.pval <- sapply(pathway.ids, function(x){
#         index <- gene.data$go_id==x
#         pathway.genes <- gene.data[index,]$ensembl_gene_id
#         index <- gene.subset$ensembl_gene_id %in% pathway.genes
#         loadings.in_pathway <- gene.subset[index,]$avg_loading
#         loadings.out_of_pathway <- gene.subset[!index,]$avg_loading
#         if (length(pathway.genes) >= 6 # Question: why is this 6 and is this necessary and does it vary by different tests?
#             && length(loadings.in_pathway)>0
#             && length(loadings.out_of_pathway)>0
#         ){ 
#             stat.test(
#                 as.numeric(loadings.in_pathway), 
#                 as.numeric(loadings.out_of_pathway), 
#                 alternative = "greater"
#                 )$p.val
#         } else {
#             10 # set as 10>1 so as to be omitted after sorting 
#         }
#     })
#     enrichment.result$go_id <- pathway.ids
#     enrichment.result$p_value <- pathway.pval
#     return(enrichment.result)
# }

# sort the enrichment result by p-value,
# then output the terms with the <num_terms> smallest p-values

enrichment.output.greater <- function(enrichment.result, gene.data, num_terms, filename, p_cutoff=0.05){
    pval <- enrichment.result$p_value.greater
    num_terms <- min(
        num_terms,
        length(which(sort(pval)<=p_cutoff))
        )
    if (num_terms==0){
        return(data.frame())
    }  
    output.data <- data.frame(
        # information included in output
        p_value.greater=array("",num_terms),
        genes.in.set.count=array("",num_terms),
        genes.out_of.set.count=array("",num_terms),
        go_id=array("",num_terms),
        name_1006=array("",num_terms),
        go_linkage_type=array("",num_terms),
        namespace_1003=array("",num_terms),
        definition_1006=array("",num_terms),
        stringsAsFactors = F
    )
    output.data$p_value.greater <- sort(pval)[1:num_terms]
    output.data$go_id <- enrichment.result$go_id[order(pval)[1:num_terms]]
    output.data$genes.in.set.count <- enrichment.result$genes.in.set.count[order(pval)[1:num_terms]]
    output.data$genes.out_of.set.count <- enrichment.result$genes.out_of.set.count[order(pval)[1:num_terms]]
    for (num_id in 1:num_terms){
        id<-output.data$go_id[num_id]
        output.data$name_1006[num_id] <- unique(gene.data[gene.data$go_id==id,]$name_1006)
        output.data$go_linkage_type[num_id] <- paste(unique(gene.data[gene.data$go_id==id,]$go_linkage_type),collapse = ",")
        output.data$namespace_1003[num_id] <- unique(gene.data[gene.data$go_id==id,]$namespace_1003)
        output.data$definition_1006[num_id] <- unique(gene.data[gene.data$go_id==id,]$definition_1006)
    }
    write.csv(output.data,file=paste(filename,"greater.csv",sep=""))
    return(output.data)
}

enrichment.output.less <- function(enrichment.result, gene.data, num_terms, filename, p_cutoff=0.05){
    pval <- enrichment.result$p_value.less
    num_terms <- min(
        num_terms,
        length(which(sort(pval)<=p_cutoff))
        )
    if (num_terms==0){
        return(data.frame())
    }  
    output.data <- data.frame(
        # information included in output
        p_value.less=array("",num_terms),
        genes.in.set.count=array("",num_terms),
        genes.out_of.set.count=array("",num_terms),
        go_id=array("",num_terms),
        name_1006=array("",num_terms),
        go_linkage_type=array("",num_terms),
        namespace_1003=array("",num_terms),
        definition_1006=array("",num_terms),
        stringsAsFactors = F
    )
    output.data$p_value.less <- sort(pval)[1:num_terms]
    output.data$go_id <- enrichment.result$go_id[order(pval)[1:num_terms]]
    output.data$genes.in.set.count <- enrichment.result$genes.in.set.count[order(pval)[1:num_terms]]
    output.data$genes.out_of.set.count <- enrichment.result$genes.out_of.set.count[order(pval)[1:num_terms]]
    for (num_id in 1:num_terms){
        id<-output.data$go_id[num_id]
        output.data$name_1006[num_id] <- unique(gene.data[gene.data$go_id==id,]$name_1006)
        output.data$go_linkage_type[num_id] <- paste(unique(gene.data[gene.data$go_id==id,]$go_linkage_type),collapse = ",")
        output.data$namespace_1003[num_id] <- unique(gene.data[gene.data$go_id==id,]$namespace_1003)
        output.data$definition_1006[num_id] <- unique(gene.data[gene.data$go_id==id,]$definition_1006)
    }
    write.csv(output.data,file=paste(filename,"less.csv",sep=""))
    return(output.data)
}

enrichment.output.two.sided <- function(enrichment.result, gene.data, num_terms, filename, p_cutoff=0.05){
    pval <- enrichment.result$p_value.two.sided
    num_terms <- min(
        num_terms,
        length(which(sort(pval)<=p_cutoff))
        )
    if (num_terms==0){
        return(data.frame())
    }  
    output.data <- data.frame(
        # information included in output
        p_value.two.sided=array("",num_terms),
        genes.in.set.count=array("",num_terms),
        genes.out_of.set.count=array("",num_terms),
        go_id=array("",num_terms),
        name_1006=array("",num_terms),
        go_linkage_type=array("",num_terms),
        namespace_1003=array("",num_terms),
        definition_1006=array("",num_terms),
        stringsAsFactors = F
    )
    output.data$p_value.two.sided <- sort(pval)[1:num_terms]
    output.data$go_id <- enrichment.result$go_id[order(pval)[1:num_terms]]
    output.data$genes.in.set.count <- enrichment.result$genes.in.set.count[order(pval)[1:num_terms]]
    output.data$genes.out_of.set.count <- enrichment.result$genes.out_of.set.count[order(pval)[1:num_terms]]
    for (num_id in 1:num_terms){
        id<-output.data$go_id[num_id]
        output.data$name_1006[num_id] <- unique(gene.data[gene.data$go_id==id,]$name_1006)
        output.data$go_linkage_type[num_id] <- paste(unique(gene.data[gene.data$go_id==id,]$go_linkage_type),collapse = ",")
        output.data$namespace_1003[num_id] <- unique(gene.data[gene.data$go_id==id,]$namespace_1003)
        output.data$definition_1006[num_id] <- unique(gene.data[gene.data$go_id==id,]$definition_1006)
    }
    write.csv(output.data,file=paste(filename,"two.sided.csv",sep=""))
    return(output.data)
}

enrichment.output <- function(enrichment.result, gene.data, num_terms, filename, p_cutoff=0.05){
        enrichment.output.greater(enrichment.result, gene.data, num_terms, filename, p_cutoff)
        enrichment.output.less(enrichment.result, gene.data, num_terms, filename, p_cutoff)
        enrichment.output.two.sided(enrichment.result, gene.data, num_terms, filename, p_cutoff)
}

enrichment.output.all <- function(enrichment.result, gene.data, filename){
    index <- (enrichment.result$p_value.larger<=1 &&
              enrichment.result$p_value.less<=1 &&
              enrichment.result$p_value.two.sided<=1)
    enrichment.result<-enrichment.result[index,]
    num_terms <- length(enrichment.result$go_id)

    output.data <- data.frame(
        # information included in output
        p_value.larger=array("",num_terms),
        p_value.less=array("",num_terms),
        p_value.two.sided=array("",num_terms),
        genes.in.set.count=array("",num_terms),
        genes.out_of.set.count=array("",num_terms),
        go_id=array("",num_terms),
        name_1006=array("",num_terms),
        go_linkage_type=array("",num_terms),
        namespace_1003=array("",num_terms),
        definition_1006=array("",num_terms),
        stringsAsFactors = F
    )

    output.data$p_value.larger <- enrichment.result$p_value.larger
    output.data$p_value.less <- enrichment.result$p_value.less
    output.data$p_value.two.sided <- enrichment.result$p_value.two.sided
    output.data$genes.in.set.count <- enrichment.result$genes.in.set.count
    output.data$genes.out_of.set.count <- enrichment.result$genes.out_of.set.count
    output.data$go_id <- enrichment.result$go_id

    for (num_id in 1:num_terms){
        id<-output.data$go_id[num_id]
        output.data$name_1006[num_id] <- unique(gene.data[gene.data$go_id==id,]$name_1006)
        output.data$go_linkage_type[num_id] <- paste(unique(gene.data[gene.data$go_id==id,]$go_linkage_type),collapse = ",")
        output.data$namespace_1003[num_id] <- unique(gene.data[gene.data$go_id==id,]$namespace_1003)
        output.data$definition_1006[num_id] <- unique(gene.data[gene.data$go_id==id,]$definition_1006)
    }
    
    write.csv(output.data,file=paste(filename,".ALL.csv",sep=""))
}