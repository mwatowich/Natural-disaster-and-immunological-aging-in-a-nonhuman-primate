### Linear modeling of whole blood RNA seq counts 

# Load normalized gene counts 
norm_gene_counts <- read.delim("path/to/dir")

# Load metadata and make design matrix 
metadata<-as.data.frame(read_xlsx('~/Downloads/PNAS_202121663_s2.xlsx', skip = 1, col_names = T))
d_matrix <- model.matrix(~ metadata$Sex + metadata$`Chronological age`, metadata$`Hurricane exposure`, metadata$`RNA quality score`)

# Load identity matrix 
Zmatrix<-read.delim("~/Documents/_UW/_SMackLab/Hurricane/Github/data/Zmatrix.csv",sep = ",",check.names = F)
rownames(Zmatrix)<-Zmatrix[,1]
Zmatrix[,1]<-NULL

# Load kinship matrix 
kinship<-read.delim("~/Documents/_UW/_SMackLab/Hurricane/Github/data/kinship.csv",sep = ",",check.names = F)
rownames(kinship)<-kinship[,1]
kinship[,1]<-NULL

# Linear modeling using EMMREML
library(EMMREML)
rna_mod_out <- data.frame(t(apply(norm_gene_counts, 1, function(y){
  emma <- emmreml(y = y,
                  X = d_matrix,
                  Z = as.matrix(Zmatrix),
                  K = as.matrix(kinship), 
                  varbetahat = T, varuhat = T, PEVuhat = T, test = T)
  p <- emma$pvalbeta
  se <- sqrt(emma$varbetahat)
  b <- emma$betahat
  return(c(b[-1,], se[-1], p[-1, "none"]))
})))

