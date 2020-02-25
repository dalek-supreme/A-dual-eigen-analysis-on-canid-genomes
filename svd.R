# Load the R packages: gdsfmt and SNPRelate
library(gdsfmt)
library(SNPRelate)

#vcf.fn <- "/Users/leo/genome-data/simp_2012.5.18_6canid_merged_chr01_cf31.vcf"
#Reformat vcf to gds, use bialletic only for this try
#Actually there are some non-bi-allelic sites, deal with this later.
#snpgdsVCF2GDS(vcf.fn, "test.gds", method="biallelic.only")

#Summary:
#snpgdsSummary("test.gds")
#Open gds file and check some parameters

genofile <- snpgdsOpen("/Users/leo/test.gds")

#get.attr.gdsn(index.gdsn(genofile,"snp.chromosome"))
#Get genotypical matrix

# find which chromosome
chrome <- read.gdsn(index.gdsn(genofile, "snp.chromosome"))
# find which sample
samples <- read.gdsn(index.gdsn(genofile, "sample.id"))
# get the CanFam 3.1 coordinates
pos <- read.gdsn(index.gdsn(genofile, "snp.position"))

ntypes<-1
chr_type<-chrome[1]
sizes<-1
for (i in 2:length(pos)){
	if (chrome[i]!=chr_type[ntypes]){
		ntypes<-ntypes+1
		chr_type[ntypes]<-chrome[i]
		sizes[ntypes]<-i
	}
}
ntypes<-ntypes+1
sizes[ntypes]<-i+1

###### 对每一条染色体做，循环
# for (chr_t in 1:(ntypes-1)){
###### 现在先用第25个，实际是31号
# !!!!! 这里的chr_t并不是染色体的编号，而是在gds文件中出现的顺序，实际编号为chrome[sizes[chr_t]]

chr_t<-25

g <- read.gdsn(index.gdsn(genofile, "genotype"), start=c(1,sizes[chr_t]), count=c(length(samples),sizes[chr_t+1]-sizes[chr_t]))
# why 998500?
#(g <- read.gdsn(index.gdsn(genofile, "genotype"), start=c(1,1), count=c(6,80)))

#################################### data cleaning #####################################
# Deal with missing genotypes, namingly the '3's in the g matrix
# first try: simply discard the columns containing '3's
# also, switch 0 and 2 to make more sense

###########################################################################################################
# is the 3 missing in the experiment or no allele? if no allele, should use -1 to make sense?  i.e. INDEL?#
###########################################################################################################

dimensions <- dim(g)
data <- array(-1, dim=dimensions)
coordinates <- array(-1,dim=dimensions[2])

j <- 1
coord <- 1
while (coord <= dimensions[2]){
  has_data <- TRUE # has data
  i <- 1
  while (i <= dimensions[1]){
    if (g[i,coord]==0){
      data[i,j] <- 2
    }
    else if (g[i,coord]==2){
      data[i,j] <- 0
    }
    else if (g[i,coord]==1){
      data[i,j] <- 1
    }
    else { # =3
      coord <- coord+1
      has_data <- FALSE # has no data (missing)
      break
    }
    i <- i+1
  }
  if (has_data) {
    coordinates[j]=coord
    j <- j+1
    coord <- coord+1
  }
}
size<-1
while (size <= dimensions[2]){
 if (coordinates[size]==-1){
  break 
 }
 size <- size+1
}
size <- size-1
coordinates <- coordinates[1:size]
data<-data[,1:size]
#valid_rate[chr_t]<-size/dimensions[2] # how much of data is good and left for further analysis

############################### end of data cleaning ####################################

############################### SVD ####################################
s <- svd(data)
#D<-diag(s$d)
#U<-s$u
#V<-s$v

## check if the svd is done properly
#A<-U %*% D %*% t(V)
#dim_valid=dim(data)
#l<-0
#for (i in 1:dim_valid[1]){
#  for (j in 1:dim_valid[2]){
#    if((A[i,j]-data[i,j]>1e-9)){
#      print("NNNNOOOOOO")
#      print(i)
#      print(j)
#      print(A[i,j]-data[i,j])
#      l<-l+1
#    }
#  }
#}
#print(l)
############################# end of SVD ################################
#print(s$d)
#}

 