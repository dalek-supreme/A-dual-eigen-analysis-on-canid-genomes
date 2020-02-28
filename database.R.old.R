# Obtain Database for GO enrichment
require(AnnotationHub)
hub <- AnnotationHub()
# hub[hub$species=='Canis familiaris'&hub$rdataclass=='OrgDb']
sql_id_canfam<-names(hub[hub$species=='Canis familiaris'&hub$rdataclass=='OrgDb'])
# Canisfam<-hub[["AH75736"]]
Canisfam<-hub[[sql_id_canfam]]

GO.data<-select(Canisfam,keys(Canisfam),c("GENENAME","GO"))
save(GO.data,file="GO.data.Rdata")

# Obtain Gene ENTRENZIDs from SNP sites
# initiallize get mart:
library("biomaRt")
ensembl = useMart("ensembl",dataset="cfamiliaris_gene_ensembl")
gene.id <- getBM(attributes=c('entrezgene_id'),
	 		filters=c('chromosomal_region'),
	 		values=paste(
                 chrome[sizes[chr_t]],
                 pos[coordinates[order(s$v[,1])[1:100]]],
                 pos[coordinates[order(s$v[,1])[1:100]]]
                 ,sep=':'),
	 		mart=ensembl)

testa <- 0
for (i in 1:10){
    testa[i] <- getBM(attributes=c('entrezgene_id'),
	     		filters=c('chromosomal_region'),
	 	    	values=paste(
                     chrome[sizes[chr_t]],
                     pos[coordinates[order(s$v[,1])[i]]],
                     pos[coordinates[order(s$v[,1])[i]]]
                     ,sep=':'),
	 		    mart=ensembl)
    print(paste(i,testa[i],sep=","))
}


#####################################################
# NEW WAY: FROM BIOMART, SUM LOADINGS BY EACH GENE! #
#####################################################
library("biomaRt")
ensembl = useMart("ensembl",dataset="cfamiliaris_gene_ensembl")
gene.data <- getBM(attributes=c("start_position","end_position","go_id","name_1006","definition_1006","go_linkage_type","namespace_1003"),
	 		filters='chromosome_name',
	 		values=chr_t,
	 		mart=ensembl)
# for chromosome=chr_t