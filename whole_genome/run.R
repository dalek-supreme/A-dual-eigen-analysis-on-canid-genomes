# pre-processing genes before enrichment

gene.loading.stats <- function(snp.svd.right, gene.data) {
    gene.isunique <- (!duplicated(gene.data$ensembl_gene_id)
                  & gene.data$go_id!="") # remove genes w/o GO annotation
    gene.unique.index <- seq_along(gene.data$ensembl_gene_id)[gene.isunique]
    gene.unique <- data.frame(
                    ensembl_gene_id = gene.data$ensembl_gene_id[gene.unique.index],
                    chromosome_name=gene.data$chromosome_name[gene.unique.index],
                    start_position = gene.data$start_position[gene.unique.index],
                    end_position = gene.data$end_position[gene.unique.index],
                    stringsAsFactors = F
                )
    
    result <- list(
        gene.unique$ensembl_gene_id, # "ensembl_gene_id"
        array(0,length(gene.unique$ensembl_gene_id)), # "loading.count"
        array(0,c(length(gene.unique$ensembl_gene_id),dim(snp.svd.right)[2])), # "loading.sum"
        array(0,c(length(gene.unique$ensembl_gene_id),dim(snp.svd.right)[2])) # "loading.avg"
    )
    names(result) <- c("ensembl_gene_id","loading.count","loading.sum","loading.avg")

    total.genes <- length(gene.unique$ensembl_gene_id)
    total.sites <- dim(snp.svd.right)[1]

    for (gene.index in seq_along(gene.unique$ensembl_gene_id)) {
        gene.loading.count <- 0
        gene.loading.sum <- array(0,dim(snp.svd.right)[2]) # do it once for all layers
        gene.loading.avg <- array(0,dim(snp.svd.right)[2])
        for (site.index in seq(dim(snp.svd.right)[1])) {
            if (snp.list$chromosome[site.index] == gene.unique$chromosome_name[gene.index]){
                if (gene.unique$start_position[gene.index] <= snp.list$position[site.index]
                && gene.unique$end_position[gene.index] >= snp.list$position[site.index]
                ) {
                    gene.loading.count <- gene.loading.count+1
                    gene.loading.sum <- gene.loading.sum + snp.svd.right[site.index,]
                    
                    #print(paste("gene:",gene.index,"/",total.genes,"=",gene.index/total.genes,sep=""))
                    #print(paste("site:",site.index,"/",total.sites,"=",site.index/total.sites,sep=""))
                    # t2<-proc.time()-t1
                    # print(t2[3])
                    # print(paste("Estimated time remaining: ",t2[3]*total.genes/gene.index-t2[3]," seconds",sep=""))
                }
            }
        }
        if (gene.loading.count != 0) {
            gene.loading.avg <- gene.loading.sum/gene.loading.count
        }
        result$loading.count[gene.index] <- gene.loading.count
        result$loading.sum[gene.index,] <- gene.loading.sum
        result$loading.avg[gene.index,] <- gene.loading.avg

        print(paste("gene:",gene.index,"/",total.genes,"=",gene.index/total.genes,sep=""))
        t2<-proc.time()-t1
        print(t2[3])
        print(paste("Estimated time remaining: ",t2[3]*total.genes/gene.index-t2[3]," seconds",sep=""))
    }
    return(result)
}

# returns genes and corresponding non-zero average loadings
gene.loadings.stats.get_layer <- function(gene.stats,layer){
    index.nonzero <- seq_along(gene.stats$loading.count)[gene.stats$loading.count!=0]
    result <- data.frame(
        ensembl_gene_id=gene.stats$ensembl_gene_id[index.nonzero],
        loading.avg=gene.stats$loading.avg[index.nonzero,layer],
        stringsAsFactors = F
    )
    return(result)
}

gene.stats.sort <- function(gene.stats.layer){
    result <- data.frame(
        ensembl_gene_id=array("",length(gene.stats.layer$ensembl_gene_id)),
        loading.avg=array(0,length(gene.stats.layer$loading.avg)),
        stringsAsFactors = F
    )
    avg <- gene.stats.layer$loading.avg
    result$loading.avg <- sort(avg)
    result$ensembl_gene_id <- 
        gene.stats.layer$ensembl_gene_id[order(avg)]
    return(result)
}

gene.subset.select <- function(gene.stats.sorted, index.range) {
    if (NA %in% gene.stats.sorted[index.range,1]){
        stop(paste("get.gene.subset: index.start or index.end out of bound! Legit range: ",range(gene.stats.sorted[,1]),sep=""))
    }
    return(gene.stats.sorted[index.range,])
}

gene.subset.pole_split <- function(gene.stats.layer,pole){
    if (pole=="positive"){
        return(gene.stats.layer[gene.stats.layer$loading.avg>0,])
    } else if (pole=="negative") {
        return(gene.stats.layer[gene.stats.layer$loading.avg<0,])
    } else {
        stop(paste(" Invalid argument \"",pole,"\": should be \"positive\" or \"negative\"!\n",sep=""))
    }
}

# enrichment analysis and file output

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
    #  t1 <- proc.time()
    pathway.ids <- unique(gene.data$go_id)
    #  pathway.total <- length(pathway.ids)
    if (length(gene.subset$ensembl_gene_id)<6) {
        enrichment.result$p_value.greater <- array(10,length(unique(gene.data$go_id)))# set as 10>1 so as to be omitted after sorting 
        enrichment.result$p_value.less <- array(10,length(unique(gene.data$go_id)))
        enrichment.result$p_value.two.sided <- array(10,length(unique(gene.data$go_id)))
        enrichment.result$go_id <- pathway.ids
        return(enrichment.result)
    }
    for (id_num in seq_along(pathway.ids)){
        index <- (gene.data$go_id==pathway.ids[id_num])
        pathway.genes <- gene.data[index,]$ensembl_gene_id
        index <- gene.subset$ensembl_gene_id %in% pathway.genes
        loadings.in_pathway <- gene.subset[index,]$loading.avg
        loadings.out_of_pathway <- gene.subset[!index,]$loading.avg
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
    #  print(paste("gene:",id_num,"/",pathway.total,"=",id_num/pathway.total,sep=""))
    #  t2 <- proc.time()-t1
    #  print(paste("Time elapsed in this layer: ",t2[3],". Estimated time remaining: ",t2[3]*pathway.total/id_num-t2[3],sep=""))
    }

    enrichment.result$go_id <- pathway.ids
    return(enrichment.result)
}

# sort the enrichment result by p-value,
# then output the terms with the <num_terms> smallest p-values

enrichment.output.greater <- function(enrichment.result, gene.data, num_terms=20, filename="wilcox", p_cutoff=0.05){
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
    write.csv(output.data,file=paste(filename,".greater.csv",sep=""))
    return(output.data)
}

enrichment.output.less <- function(enrichment.result, gene.data, num_terms=20, filename="wilcox", p_cutoff=0.05){
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
    write.csv(output.data,file=paste(filename,".less.csv",sep=""))
    return(output.data)
}

enrichment.output.two.sided <- function(enrichment.result, gene.data, num_terms=20, filename="wilcox", p_cutoff=0.05){
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
    write.csv(output.data,file=paste(filename,".two.sided.csv",sep=""))
    return(output.data)
}

enrichment.output <- function(enrichment.result, gene.data, num_terms=20, filename="wilcox", p_cutoff=0.05){
        enrichment.output.greater(enrichment.result, gene.data, num_terms, filename, p_cutoff)
        enrichment.output.less(enrichment.result, gene.data, num_terms, filename, p_cutoff)
        enrichment.output.two.sided(enrichment.result, gene.data, num_terms, filename, p_cutoff)
}



load("gene.db.Rdata")
load("gene.stats.Rdata")

num_layers <- 10
layers <- seq(num_layers)
num_terms <- 50

t1 <- proc.time()
print("Starting gene enrichment analysis...")
t5<-proc.time()
for (layer in layers)
{
    t3 <- proc.time()
    gene.stats.layer <- gene.loadings.stats.get_layer(gene.stats,layer)
    gene.stats.sorted <-  gene.stats.sort(gene.stats.layer)
    gene.stats.positive <- gene.subset.pole_split(gene.stats.sorted,pole="positive")
    gene.stats.negative <- gene.subset.pole_split(gene.stats.sorted,pole="negative")

    result.positive <- enrichment.wilcoxon(gene.stats.positive,gene.db) # about 7 minutes
    result.negative <- enrichment.wilcoxon(gene.stats.negative,gene.db)
    setwd(paste("~/enrichment/",layer,sep=""))
    save(result.positive,result.negative,file="result.Rdata")
    setwd(paste("~/enrichment/",layer,"/positive",sep=""))
    enrichment.output(result.positive, gene.db, num_terms, filename="wilcox", p_cutoff=0.05)
    setwd(paste("~/enrichment/",layer,"/negative",sep=""))
    enrichment.output(result.negative, gene.db, num_terms, filename="wilcox", p_cutoff=0.05)
    #enrichment.output.all(result, gene.data, filename)
    t4 <- proc.time()-t3
    print(paste("Layer ",layer," finished. Time elapsed: ",t4[3],sep=""))
    t6<-proc.time()-t5
    print(paste("Estimated time remaining: ",t6[3]*(length(layers)/layer-1)/60,"minutes.",sep=""))
}

t7 <- proc.time()-t0
print(paste("Analysis complete. Total time: ",t7[3],sep=""))

#23571