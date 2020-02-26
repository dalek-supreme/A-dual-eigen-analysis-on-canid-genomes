# Obtain database of gene, including:
# Gene start/end position,
# GO accession code,
# GO term name,
# GO term description,
# and GO evidence.
# currently just for one chromosome.
library("biomaRt")
ensembl = useMart("ensembl",dataset="cfamiliaris_gene_ensembl")
gene.data <- getBM(attributes=c("start_position","end_position","go_id","name_1006","definition_1006","go_linkage_type","namespace_1003"),
	 		filters='chromosome_name',
	 		values=chr_t,
	 		mart=ensembl)
# for chromosome=chr_t