# library(AnnotationHub)
# hub <- AnnotationHub()
# # hub[hub$species=='Canis familiaris'&hub$rdataclass=='OrgDb']
# sql_id_canfam<-names(hub[hub$species=='Canis familiaris'&hub$rdataclass=='OrgDb'])
# # Canisfam<-hub[["AH75736"]]
# Canisfam<-hub[[sql_id_canfam]]

# egid <- head(keys(Canisfam, "ENTREZID"))
# select(Canisfam, egid, c("SYMBOL", "PMID"), "ENTREZID")

library(clusterProfiler)
KEGG.db <- download_KEGG(species='cfa')
save(KEGG.db,'KEGG.db.Rdata')