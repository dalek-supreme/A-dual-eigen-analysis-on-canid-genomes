#install packages:
# SNPRelate
# biomaRt
# gdsfmt
# clusterProfiler
# AnnotationHub
# rtracklayer

# if (!require(SNPRelate)) install.packages('SNPRelate')
# library(SNPRelate)

# if (!require(biomaRt)) install.packages('biomaRt')
# library(biomaRt)

# if (!require(gdsfmt)) install.packages('gdsfmt')
# library(gdsfmt)

# if (!require(clusterProfiler)) install.packages('clusterProfiler')
# library(clusterProfiler)

# if (!require(AnnotationHub)) install.packages('AnnotationHub')
# library(AnnotationHub)

# if (!require(rtracklayer)) install.packages('rtracklayer')
# library(rtracklayer)

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("SNPRelate")
BiocManager::install("biomaRt")
BiocManager::install("gdsfmt")

#BiocManager::install("pcaMethods")
#install.packages("matrixStats")

install.packages("RobRSVD")
install.packages("rsvd")
install.packages("ape")
install.packages("tidytree")
install.packages("phytools")

install.packages("ggplot2")
install.packages("ggmap")
install.packages("sp")
install.packages("maptools")
install.packages("maps")
install.packages("wesanderson")