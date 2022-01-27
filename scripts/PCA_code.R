# Perform PCA of global gene expression 

# Remove partial residuals of RQN as we know that this technical variance will affect expression
# Load linear model results and subtract out partial residuals of RQN (normalized_GE - Brqn*rqn)
resid_GE<-norm_gene_counts-(as.matrix(rna_mod_out$beta_RQN)%*%t(as.matrix(metadata$RQN)))

# Perform PCA and model 
counts_pca = prcomp(cor(resid_GE), center = T, scale = T)
summary(counts_pca)[["importance"]] %>% 
  as.data.frame() %>% 
  dplyr::select(c(PC1, PC2, PC3))  
counts_pca <- as.data.frame(counts_pca$x) %>% 
  rownames_to_column(var = "LID")
counts_pc_meta <- left_join(counts_pca, sample_info, by = "LID")
summary(lm(PC1 ~ sex + hurricane + FA_RQN + age_at_blood, counts_pc_meta))
