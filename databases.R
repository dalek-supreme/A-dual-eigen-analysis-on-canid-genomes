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
                 pos[coordinates],
                 pos[coordinates]
                 ,sep=':'),
	 		mart=ensembl)