### Normalize raw gene counts

# Load read counts 
gene_counts <- read.delim('path/to/gene_counts.csv',sep = '\t')
rownames(gene_counts)<-gene_counts$X
gene_counts$X<-NULL

# Normalize to counts per million
cpm = as.data.frame(apply(gene_counts,2,function(x){x/sum(x)*1e6}))

# Filter out lowly expressed genes (plot to find cut-off)
cpm <- cpm[apply(cpm, 1, function(x){median(x) > 3.277}),]

# Get raw counts for only have the genes that passed threshold 
filtered_gene_counts <- gene_counts[rownames(cpm),]

# Voom normalize filtered counts 
library(limma)
norm_gene_counts <- voom(counts = filtered_gene_counts)$E
