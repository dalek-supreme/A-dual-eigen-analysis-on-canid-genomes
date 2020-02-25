require(AnnotationHub)
hub <- AnnotationHub()
# hub[hub$species=='Canis familiaris'&hub$rdataclass=='OrgDb']
sql_id_canfam<-names(hub[hub$species=='Canis familiaris'&hub$rdataclass=='OrgDb'])
# Canisfam<-hub[["AH75736"]]
Canisfam<-hub[[sql_id_canfam]]


# translate ID to Ontology terms
#library(clusterProfiler)
#bitr(keys(Canisfam)[1], 'ENTREZID', c("REFSEQ", "GO", "ONTOLOGY"), Canisfam)

#save(Canisfam,file="Canisfam.Rdata")