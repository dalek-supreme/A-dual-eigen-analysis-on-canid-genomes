# "DualEigen_utils.R"
principle.loadings <- function(vec,percentile=0.9,number,UseNumber=F) {
    if (!UseNumber){
        result <- vector()
        sval <- quantile(vec,1-percentile)
        lval <- quantile(vec,percentile)
        for (val in vec){
            if (val<sval | val>lval) {
                result<-append(result,val)
            }
        }
        return(sort(result))
    }
    else {
        sorted <- sort(vec)
        result <- vector()
        for (i in seq(min(number/2,length(vec)))) {
            result <- append(result,sorted[i])
        }
        for (i in length(vec)-seq(min(number/2,length(vec)))) {
            result <- append(result,sorted[i])
        }
        return(sort(result))
    }
}

non.principle.loadings <- function(vec,percentile=0.9,number,UseNumber=F) {
    if (!UseNumber){
        result <- vector()
        sval <- quantile(vec,1-percentile)
        lval <- quantile(vec,percentile)
        for (val in vec){
            if (val>sval & val<lval) {
                result<-append(result,val)
            }
        }
        return(sort(result))
    }
    else {
        sorted <- sort(vec)
        result <- sorted[seq(number/2,length(vec)-number/2)]
        return(result)
    }
}

filter.loadings <- function(vec,percentile=0.9,norm.cutoff=0.9,number=20,threshold,method=c('2','percentile','number','threshold')) {
    result <- array(0,length(vec))
    if (method == 'percentile'){
        sval <- quantile(vec,1-percentile)
        lval <- quantile(vec,percentile)
        for (i in seq_along(vec)){
            if (vec[i]<sval | vec[i]>lval) {
                result[i]<-vec[i]
            }
        }
    }
    else if (method == 'number'){
        threshold <- sort(abs(vec))[length(vec)-number+1]
        for (i in seq_along(vec)){
            if (abs(vec[i])>=threshold) {
                result[i]<-vec[i]
            }
        }
    }
    else if (method == '2'){
        sorted <- sort(vec^2,decreasing=T)
        sum <- 0
        vec.norm.sqr <- sum(vec^2)
        for (value in sorted){
            sum <- sum + value
            if (sum >= vec.norm.sqr*norm.cutoff^2) break;
        }
        threshold <- sqrt(value)
        for (i in seq_along(vec)){
            if (abs(vec[i])>=threshold) {
                result[i]<-vec[i]
            }
        }
    }
    else if (method=='threshold'){
        for (i in seq_along(vec)){
            if (abs(vec[i])>=threshold) {
                result[i]<-vec[i]
            }
        } 
    }
    else stop('invalid argument.\n')
    return(result)
}

nonzero.pos <- function(vec) {
    result <- vector()
    for (i in seq_along(vec)){
        if (vec[i]!=0) result <- append(result,i)
    }
    return(result)
}

random.pos <- function(total.length,size){
    return(sample(seq(total.length),size=size,replace=F))
}


unit.vec.angle <- function(v1,v2,unit=c('rad','deg')){
    if (unit == 'rad') return(acos(sum(v1*v2)))
    if (unit == 'deg') return(acos(sum(v1*v2))/pi*180)
    else stop('Wrong argument for unit: use \'rad\' or \'deg\'\n')
}


# "utils-2.R"
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

enrichment.wilcoxon <- function(gene.subset,gene.data, pathway.namespace){
    gene.subset <- gene.subset[gene.subset$ensembl_gene_id %in% unique(gene.data[gene.data$namespace_1003==pathway.namespace,]$ensembl_gene_id),]
    pathway.ids <- unique(gene.data[gene.data$namespace_1003==pathway.namespace,]$go_id)
    result.length <- length(pathway.ids)
    enrichment.result <- data.frame(
        go_id = array("",result.length),
        p_value.greater = array("",result.length),
        p_value.less = array("",result.length),
        p_value.two.sided = array("",result.length),
        genes.in.set.count = array("",result.length),
        genes.out_of.set.count = array("",result.length),
        stringsAsFactors = F
    )
    #  t1 <- proc.time()
    pathway.ids <- unique(gene.data[gene.data$namespace_1003==pathway.namespace,]$go_id)
    #  pathway.total <- length(pathway.ids)
    if (length(gene.subset$ensembl_gene_id)<6) {
        enrichment.result$p_value.greater <- array(10,result.length)# set as 10>1 so as to be omitted after sorting 
        enrichment.result$p_value.less <- array(10,result.length)
        enrichment.result$p_value.two.sided <- array(10,result.length)
        enrichment.result$go_id <- pathway.ids
        return(enrichment.result)
    }
    for (id_num in seq_along(pathway.ids)){
        index <- (gene.data$go_id==pathway.ids[id_num])
        pathway.genes <- gene.data[index,]$ensembl_gene_id
        index <- gene.subset$ensembl_gene_id %in% pathway.genes
        loadings.in_pathway <- gene.subset[index,]$loading.avg
        loadings.out_of_pathway <- gene.subset[!index,]$loading.avg # need to change. need to intersect with G?
        genes.in.set.count <- length(loadings.in_pathway)
        genes.out_of.set.count <- length(loadings.out_of_pathway)

        enrichment.result$genes.in.set.count[id_num] <- genes.in.set.count
        enrichment.result$genes.out_of.set.count[id_num] <- genes.out_of.set.count

        if (length(pathway.genes) >= 6 # Question: why is this 6 and is this necessary and does it vary by different tests?
            && genes.in.set.count>=6
            && genes.out_of.set.count>=6
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
    #   print(paste("gene:",id_num,"/",pathway.total,"=",id_num/pathway.total,sep=""))
    #   t2 <- proc.time()-t1
    #   print(paste("Time elapsed in this layer: ",t2[3],". Estimated time remaining: ",t2[3]*pathway.total/id_num-t2[3],sep=""))
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


gene.loading.stats.1_layer_only <- function(snp.svd.v.layer,gene.data) {
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
    result <- data.frame(
        ensembl_gene_id=gene.unique$ensembl_gene_id, # "ensembl_gene_id"
        loading.count=array(0,length(gene.unique$ensembl_gene_id)), # "loading.count"
        loading.sum=array(0,length(gene.unique$ensembl_gene_id)), # "loading.sum"
        loading.avg=array(0,length(gene.unique$ensembl_gene_id)), # "loading.avg"
        stringsAsFactors = F
    )
    total.genes <- length(gene.unique$ensembl_gene_id)
    total.sites <- length(snp.svd.v.layer)

    for (gene.index in seq_along(gene.unique$ensembl_gene_id)) {
        gene.loading.count <- 0
        gene.loading.sum <- 0
        gene.loading.avg <- 0
        for (site.index in seq(length(snp.svd.v.layer))) {
            if (snp.list$chromosome[site.index] == gene.unique$chromosome_name[gene.index]){
                if (gene.unique$start_position[gene.index] <= snp.list$position[site.index]
                && gene.unique$end_position[gene.index] >= snp.list$position[site.index]
                ) {
                    gene.loading.count <- gene.loading.count+1
                    gene.loading.sum <- gene.loading.sum + (snp.svd.v.layer[site.index])^2 # changed to square. see if any different..
                    
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
        result$loading.sum[gene.index] <- gene.loading.sum
        result$loading.avg[gene.index] <- gene.loading.avg

        # print(paste("gene:",gene.index,"/",total.genes,"=",gene.index/total.genes,sep=""))
        # t2<-proc.time()-t1
        # print(t2[3])
        # print(paste("Estimated time remaining: ",t2[3]*total.genes/gene.index-t2[3]," seconds",sep=""))
    }
    result <- result[result$loading.count!=0,]
    return(result)
}

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
# source('utils-2.R')
# source('DualEigen_utils.R')
#load('snp.svd.Rdata')
#load('snp.genotype.Rdata')

#sink("WILCOXON_enrichment.log")

# Wilcoxon Rank Sum Enrichment
num_layers <- 10
layers <- 2:10
num_terms <- 500

system('cd ~')
system('mkdir enrichment_GO_all_square')
# system(paste("NUM_LAYERS=",num_layers,sep=""))
# system("bash create_result_dir.sh")

t0 <- proc.time()
t1 <- proc.time()

# gene.stats <- gene.loading.stats(snp.svd$v,gene.db)
# ########## stuck here ##########
# save(gene.stats,file='gene.stats2.Rdata')

# ## Get main positions:
# u.cutoff <- 0.8 -> v.cutoff
# main.pos.x <- list()
# main.pos.y <- list()
# for (layer in seq(num_layers)) {
#     u <- filter.loadings(snp.svd$u[,layer],method='2',norm.cutoff=u.cutoff)
#     v <- filter.loadings(snp.svd$v[,layer],method='2',norm.cutoff=v.cutoff)
#     main.pos.x[[layer]] <- nonzero.pos(u)
#     main.pos.y[[layer]] <- nonzero.pos(v)
# }

# for (layer in seq(num_layers)) {
#     print(paste(layer,length(main.pos.x[[layer]]),length(main.pos.y[[layer]])))
# }

# print('main positions found.')

# load("gene.db.Rdata")

t5<-proc.time()
for (layer in layers)
{
    t3 <- proc.time()
    gene.stats <- gene.loading.stats.1_layer_only(snp.svd$v[,layer],gene.db)
    
    #gene.stats.layer <- gene.loadings.stats.get_layer(gene.stats)
    #gene.stats.sorted <-  gene.stats.sort(gene.stats.layer)
    # gene.stats.positive <- gene.subset.pole_split(gene.stats.sorted,pole="positive")
    # gene.stats.negative <- gene.subset.pole_split(gene.stats.sorted,pole="negative")
    #gene.stats.main <- gene.subset.select(gene.stats.layer, main.pos.y[[layer]])
    result.molecular_function <- enrichment.wilcoxon(gene.stats,gene.db,'molecular_function') # about 7 minutes
    result.biological_process <- enrichment.wilcoxon(gene.stats,gene.db,'biological_process') # about 7 minutes
    result.cellular_component <- enrichment.wilcoxon(gene.stats,gene.db,'cellular_component') # about 7 minutes
    
    setwd("~/enrichment_GO_all_square/")
    save(gene.stats,result.molecular_function,result.biological_process,result.cellular_component,file=paste('result-',layer,'.Rdata',sep=''))
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

