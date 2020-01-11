# Problem: now some layers have no good enrichment results.
# Question: how large is a good enrichment set? -- see paper.
# Solution: try a dynamic repeat until !is.na strategy, adjust gene set size

# also, try with a KEGG and a Wilcoxon



######################################################################
#                   enrichment and write to file                     #
######################################################################
# initiallize get mart:
library("biomaRt")
library('clusterProfiler')
ensembl = useMart("ensembl",dataset="cfamiliaris_gene_ensembl")

# set variables
layer<-1
gpos_size<-100

pvalC<-0.2
qvalC<-0.2

pval_range<-10

output_enrichment<-array('',dim=c(pval_range+pval_range+1+pval_range+1+pval_range+1,9))

for(layer in 8:10){
	g_orders <- order(s$v[,layer])
	len<-length(s$v[,layer])
	output_enrichment<-array('',dim=c(pval_range+pval_range+1+pval_range+1+pval_range+1,9))
	for(j in 1:pval_range){
		# LARGE SIDE	
		g_pos_l<-(len-gpos_size+1):(len-1+1)

		(g_ncbi_large <- getBM(attributes=c('entrezgene_id'),
	 		filters=c('chromosomal_region'),
	 		values=paste(chrome[sizes[chr_t]],pos[coordinates[g_orders[g_pos_l]]],pos[coordinates[g_orders[g_pos_l]]],sep=':'),
	 		mart=ensembl))
		
		print(paste("layer=",layer,", j=",j,", NCBI IDs at large side: ",g_ncbi_large,sep=""))
		
		go_large_BP <- enrichGO(as.matrix(g_ncbi_large), OrgDb=Canisfam, pAdjustMethod = "BH", pvalueCutoff=pvalC, qvalueCutoff=qvalC, ont='BP')
		go_large_CC <- enrichGO(as.matrix(g_ncbi_large), OrgDb=Canisfam, pAdjustMethod = "BH", pvalueCutoff=pvalC, qvalueCutoff=qvalC, ont='CC')
		go_large_MF <- enrichGO(as.matrix(g_ncbi_large), OrgDb=Canisfam, pAdjustMethod = "BH", pvalueCutoff=pvalC, qvalueCutoff=qvalC, ont='MF')
		
		print(go_large_BP$Description)
		print(go_large_CC$Description)
		print(go_large_MF$Description)

		
		# SMALL SIDE
		g_pos_s<-1:gpos_size
		
		(g_ncbi_small <- getBM(attributes=c('entrezgene_id'),
	 		filters=c('chromosomal_region'),
	 		values=paste(chrome[sizes[chr_t]],pos[coordinates[g_orders[g_pos_s]]],pos[coordinates[g_orders[g_pos_s]]],sep=':'),
	 		mart=ensembl))
	 	
	 	print(paste("layer=",layer,", j=",j,", NCBI IDs at small side: ",g_ncbi_small,sep=""))
	 	
	 	go_small_BP <- enrichGO(as.matrix(g_ncbi_small), OrgDb=Canisfam, pAdjustMethod = "BH", pvalueCutoff=pvalC, qvalueCutoff=qvalC, ont='BP')
		go_small_CC <- enrichGO(as.matrix(g_ncbi_small), OrgDb=Canisfam, pAdjustMethod = "BH", pvalueCutoff=pvalC, qvalueCutoff=qvalC, ont='CC')
		go_small_MF <- enrichGO(as.matrix(g_ncbi_small), OrgDb=Canisfam, pAdjustMethod = "BH", pvalueCutoff=pvalC, qvalueCutoff=qvalC, ont='MF')

		print(go_small_BP$Description)
		print(go_small_CC$Description)
		print(go_small_MF$Description)
		
		# OUTPUT
		output_enrichment[j,1]<-go_large_BP$Description[1:pval_range][j]
		output_enrichment[j,2]<-go_large_BP$pvalue[1:pval_range][j]
		output_enrichment[j,3]<-go_large_BP$qvalue[1:pval_range][j]

		output_enrichment[j,4]<-go_large_CC$Description[1:pval_range][j]
		output_enrichment[j,5]<-go_large_CC$pvalue[1:pval_range][j]
		output_enrichment[j,6]<-go_large_CC$qvalue[1:pval_range][j]
		
		output_enrichment[j,7]<-go_large_MF$Description[1:pval_range][j]
		output_enrichment[j,8]<-go_large_MF$pvalue[1:pval_range][j]
		output_enrichment[j,9]<-go_large_MF$qvalue[1:pval_range][j]
		
		output_enrichment[j+pval_range+1,1]<-go_small_BP$Description[1:pval_range][j]
		output_enrichment[j+pval_range+1,2]<-go_small_BP$pvalue[1:pval_range][j]
		output_enrichment[j+pval_range+1,3]<-go_small_BP$qvalue[1:pval_range][j]
		
		output_enrichment[j+pval_range+1,4]<-go_small_CC$Description[1:pval_range][j]
		output_enrichment[j+pval_range+1,5]<-go_small_CC$pvalue[1:pval_range][j]
		output_enrichment[j+pval_range+1,6]<-go_small_CC$qvalue[1:pval_range][j]
		
		output_enrichment[j+pval_range+1,7]<-go_small_MF$Description[1:pval_range][j]
		output_enrichment[j+pval_range+1,8]<-go_small_MF$pvalue[1:pval_range][j]
		output_enrichment[j+pval_range+1,9]<-go_small_MF$qvalue[1:pval_range][j]

	}
	demog_sorted<-array(-1,dim=c(127,9))
	u_order_layer<-order(s$u[,layer])
	# read demographic info to generate a report
	demog<-as.matrix(read.csv('/Users/leo/genome-data/demographic.csv',header=TRUE))
	for (i in 1:length(samples)){
		demog_sorted[i,1]<-demog[which(demog[,2]==samples[u_order_layer[i]]),1]
		demog_sorted[i,2]<-demog[which(demog[,2]==samples[u_order_layer[i]]),2]
		demog_sorted[i,3]<-demog[which(demog[,2]==samples[u_order_layer[i]]),3]
		demog_sorted[i,4]<-demog[which(demog[,2]==samples[u_order_layer[i]]),4]
		demog_sorted[i,5]<-demog[which(demog[,2]==samples[u_order_layer[i]]),5]
		demog_sorted[i,6]<-demog[which(demog[,2]==samples[u_order_layer[i]]),6]
		demog_sorted[i,7]<-demog[which(demog[,2]==samples[u_order_layer[i]]),7]
		demog_sorted[i,8]<-demog[which(demog[,2]==samples[u_order_layer[i]]),8]
		demog_sorted[i,9]<-demog[which(demog[,2]==samples[u_order_layer[i]]),9]
	}
	for(j in 1:pval_range){
		output_enrichment[j+pval_range+1+pval_range+1,1]<-demog_sorted[j,1]
		output_enrichment[j+pval_range+1+pval_range+1,2]<-demog_sorted[j,2]
		output_enrichment[j+pval_range+1+pval_range+1,3]<-demog_sorted[j,3]
		output_enrichment[j+pval_range+1+pval_range+1,4]<-demog_sorted[j,4]
		
		output_enrichment[j+pval_range+1+pval_range+1+pval_range+1,1]<-demog_sorted[128-j,1]
		output_enrichment[j+pval_range+1+pval_range+1+pval_range+1,2]<-demog_sorted[128-j,2]
		output_enrichment[j+pval_range+1+pval_range+1+pval_range+1,3]<-demog_sorted[128-j,3]
		output_enrichment[j+pval_range+1+pval_range+1+pval_range+1,4]<-demog_sorted[128-j,4]
	}
	
	write.csv(as.data.frame(output_enrichment),file=paste('/Users/leo/Documents/毕设/output_10/',layer,'_layer_comparison_chrome_',chrome[sizes[chr_t]],'.csv',sep=''))
	
}


