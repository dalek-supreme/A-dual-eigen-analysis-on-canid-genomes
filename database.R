# Obtain database of gene, including:
# Gene start/end position,
# GO accession code,
# GO term name,
# GO term description,
# and GO evidence.
# currently just for one chromosome.
library("biomaRt")
ensembl = useMart("ensembl",dataset="cfamiliaris_gene_ensembl")
gene.db <- getBM(
			attributes=c("ensembl_gene_id","entrezgene_id",
						 "chromosome_name","start_position","end_position",
						 "go_id","name_1006","definition_1006","go_linkage_type","namespace_1003"),
	 		mart=ensembl,
			quote = "")

#setwd("") # set to the correct working dir
save(gene.db,file="gene.db.Rdata")
# for chromosome=chr_t
cat("Annotation data successfully acquired.\n")
