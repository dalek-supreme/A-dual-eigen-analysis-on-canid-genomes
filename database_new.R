library(GO.db)
library(KEGGREST)

# KEGG
gene2pathway.genes <- keggLink('cfa','pathway')
gene2pathway.pathways <- keggLink('pathway','cfa')
gene2pathway.ncbi_ids <- keggConv("ncbi-geneid",gene2pathway.genes)

gene2pathway.pathways <- gene2pathway.pathways[names(gene2pathway.ncbi_ids)]

KEGGPATHID2EXTID <- data.frame(
    from=gene2pathway.pathways,
    to=as.numeric(gsub('ncbi-geneid:','',gene2pathway.ncbi_ids)),
    string.as.factor=F
)
rownames(KEGGPATHID2EXTID) <- NULL

pathway.name <- keggList('pathway','cfa')
KEGGPATHID2NAME <- data.frame(
     from=names(pathway.name),
     to=pathway.name,
     string.as.factor=F
    )
rownames(KEGGPATHID2NAME) <- NULL

kegg.db <- list(KEGGPATHID2EXTID,KEGGPATHID2NAME)
names(kegg.db) <- c('KEGGPATHID2EXTID','KEGGPATHID2NAME')

save(kegg.db,file='kegg.db.Rdata')
