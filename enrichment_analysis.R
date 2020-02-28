# enrichment analysis on loadings

setwd("~") # set to the correct working dir

source("database.R")
#source("utils.R")
#load("gene.data.Rdata")

# Wilcoxon Rank Sum Enrichment
layers <- 1:10
num_terms <- 50

# largest 100 genes

setwd("~/results_wilcox/large100")
filename <- paste("output_wilcox_chr=",chr_t,sep="")
print("Starting Enrichment..")
t1<-Sys.time()
round_count<-1
for (layer in layers){
    t2<-Sys.time()

    gene.stats <- get.gene.snp.counts_and_sums(layer,s$v,gene.data)
    gene.data <- get.gene.data.remove0s(gene.stats,gene.data)
    gene.stats.sorted <-  get.gene.stats.sorted(gene.stats)
    gene.subset <- get.gene.subset(gene.stats.sorted,index.range=1:100)
    result <- enrichment.wilcoxon(gene.subset,gene.data)
    filename <- paste("output_wilcox_chr=",chr_t,"_layer=",layer,sep="")
    output <- enrichment.output(result, gene.data, num_terms, filename, p_cutoff=0.05)

    t3<-Sys.time()
    print(paste("Enrichment on layer ",layer," completed. Time spent in this layer: ",t3-t2," min. Total time elapsed: ",t3-t1," min.",sep=""))
    print(paste(round_count/length(layers)*100,"% completed.",sep=""))
    round_count<-round_count+1
}
t4<-Sys.time()
print(paste("Enrichment completed. Total Time elapsed: ",t4-t1," min",sep=""))

# smallest 100 genes

setwd("~/results_wilcox/small100")
filename <- paste("output_wilcox_chr=",chr_t,sep="")
print("Starting Enrichment..")
#t1<-Sys.time()
round_count<-1
for (layer in layers){
    t2<-Sys.time()

    gene.stats <- get.gene.snp.counts_and_sums(layer,s$v,gene.data)
    gene.data <- get.gene.data.remove0s(gene.stats,gene.data)
    gene.stats.sorted <-  get.gene.stats.sorted(gene.stats)
    gene.subset <- get.gene.subset(gene.stats.sorted,index.range=(dim(gene.stats.sorted)[1]-100:dim(gene.stats.sorted)[1]))
    result <- enrichment.wilcoxon(gene.subset,gene.data)
    filename <- paste("output_wilcox_chr=",chr_t,"_layer=",layer,sep="")
    output <- enrichment.output(result, gene.data, num_terms, filename, p_cutoff=0.05)

    t3<-Sys.time()
    print(paste("Enrichment on layer ",layer," completed. Time spent in this layer: ",t3-t2," min. Total time elapsed: ",t3-t1," min.",sep=""))
    print(paste(round_count/length(layers)*100,"% completed.",sep=""))
    round_count<-round_count+1
}
t4<-Sys.time()
print(paste("Enrichment completed. Total Time elapsed: ",t4-t1," min",sep=""))


# largest 50 + smallest 50
# smallest 100 genes

setwd("~/results_wilcox/large100+small100")
filename <- paste("output_wilcox_chr=",chr_t,sep="")
print("Starting Enrichment..")
#t1<-Sys.time()
round_count<-1
for (layer in layers){
    t2<-Sys.time()

    gene.stats <- get.gene.snp.counts_and_sums(layer,s$v,gene.data)
    gene.data <- get.gene.data.remove0s(gene.stats,gene.data)
    gene.stats.sorted <-  get.gene.stats.sorted(gene.stats)
    gene.subset <- get.gene.subset(gene.stats.sorted,index.range=(c((1:50),dim(gene.stats.sorted)[1]-50:dim(gene.stats.sorted)[1])))
    result <- enrichment.wilcoxon(gene.subset,gene.data)
    filename <- paste("output_wilcox_chr=",chr_t,"_layer=",layer,sep="")
    output <- enrichment.output(result, gene.data, num_terms, filename, p_cutoff=0.05)

    t3<-Sys.time()
    print(paste("Enrichment on layer ",layer," completed. Time spent in this layer: ",t3-t2," min. Total time elapsed: ",t3-t1," min.",sep=""))
    print(paste(round_count/length(layers)*100,"% completed.",sep=""))
    round_count<-round_count+1
}
t4<-Sys.time()
print(paste("Enrichment completed. Total Time elapsed: ",t4-t1," min",sep=""))

#do_enrich <- function(...,filename,range)
