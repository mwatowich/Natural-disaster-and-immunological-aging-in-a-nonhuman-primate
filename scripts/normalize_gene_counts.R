### Normalize raw gene counts

# Load read counts 
gene_counts <- read.delim('path/to/gene_counts.csv',
            sep = ",")
rownames(gene_counts)<-gene_counts$X
gene_counts$X<-NULL

# Normalize to counts per million
cpm = as.data.frame(apply(gene_counts,2,function(x){x/sum(x)*1e6}))

# Plot sum of reads per gene to assess where to split bimodal counts distribution 
ggplot(cpm, aes(x = log2(apply(cpm, 1, sum) + 1))) + 
  geom_histogram(bins = 100, fill = "grey", color = "black") + 
  labs(x = "Log 2 sum of counts per gene (cpm)", y = "Frequency") + 
  geom_vline(xintercept = 9.7, color = "blue", lty = 2, size=1)

# Filter out lowly expressed genes 
cpm <- cpm[apply(cpm, 1, function(x){median(x) > 3.277}),]

# Get raw counts for only have the genes that passed threshold 
filtered_gene_counts <- gene_counts[rownames(cpm),]

# Voom normalize filtered counts 
library(limma)
norm_gene_counts <- voom(counts = filtered_gene_counts)$E
