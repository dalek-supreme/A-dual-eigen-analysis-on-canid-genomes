require(AnnotationHub)
hub <- AnnotationHub()
# hub[hub$species=='Canis familiaris'&hub$rdataclass=='OrgDb']
sql_id_canfam<-names(hub[hub$species=='Canis familiaris'&hub$rdataclass=='OrgDb'])
# Canisfam<-hub[["AH75736"]]
Canisfam<-hub[[sql_id_canfam]]

GO.data<-select(Canisfam,keys(Canisfam),c("GO"))
save(GO.data,file="GO.data.Rdata")