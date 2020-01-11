#install packages:
# SNPRelate
# biomaRt
# gdsfmt
# clusterProfiler
# AnnotationHub
# rtracklayer

if (!require(SNPRelate)) install.packages('SNPRelate')
library(SNPRelate)

if (!require(biomaRt)) install.packages('biomaRt')
library(biomaRt)

if (!require(gdsfmt)) install.packages('gdsfmt')
library(gdsfmt)

if (!require(clusterProfiler)) install.packages('clusterProfiler')
library(clusterProfiler)

if (!require(AnnotationHub)) install.packages('AnnotationHub')
library(AnnotationHub)

if (!require(rtracklayer)) install.packages('rtracklayer')
library(rtracklayer)
