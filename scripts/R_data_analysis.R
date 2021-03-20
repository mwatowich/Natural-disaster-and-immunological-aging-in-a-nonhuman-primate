# Load libraries 
library(factoextra)
library(readxl)
library(doParallel)
library(limma)
library(corrplot)
library(topGO)
library(tidyverse)
library(biomaRt)
library(EMMREML)
library(qvalue)
library(wesanderson)
#library(devtools)
library(heatmap3)
library(ggridges)
library(ggrepel)
library(magick)
#library(RColorBrewer)
library(glmnet)

# Make an RQN-controlled residuals matrix from EMMA model results - this needs to change with new EMMA model 
resid_GE <- norm_gene_counts - 
  as.matrix(emma_pre_post$beta_FA_RQN) %*% t(as.matrix(sample_info$FA_RQN)) # resid_GE = norm_GE - Brqn*rqn - Bsex*sex (sex if want)

# Set figure theme 
tm3_theme <- theme_classic() + 
  theme(legend.position = "right",
        axis.title.x = element_text(family="Helvetica", color = "black", size = 20, face = "bold"), 
        axis.text.x = element_text(family="Helvetica", color = "black", size = 18),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_text(family="Helvetica", color = "black", size = 20, face = "bold", angle = 90, vjust = 0.5),
        axis.text.y = element_text(family="Helvetica", color = "black", size = 18), 
        legend.title = element_text(family="Helvetica", color = "black", size = 20, face = "bold"), 
        legend.text = element_text(family = "Helvetica", color = "black", size = 20, face = "bold"))

# Set colors 
color_not_sig <- wes_palette("IsleofDogs2")[1]#wes_palette("Darjeeling2", n = 1)
color_hurr <- wes_palette("Darjeeling2", n = 5)[2]
color_age <- wes_palette("Darjeeling1", n = 5)[1] #wes_palette("Cavalcanti1")[5] 
color_sex <- wes_palette("Cavalcanti1")[1]
color_both <- "purple4" # both = age+hurricane 
#color_age_sex <- wes_palette("Darjeeling1", n = 5)[4] 
#color_hurr_sex <- wes_palette("Darjeeling1", n = 5)[2]
color_pre <- "darkgreen"
color_post <- "sienna"
#


# Median age is younger after hurricane -----------
cor.test(sample_info$age_at_blood, as.numeric(sample_info$hurricane)) # neg. corr b/w age and hurr. 
tapply(sample_info$age_at_blood, sample_info$hurricane, median)

# PCA of GE: what is driving the variance in GE? ----------------------
counts_pca = prcomp(cor(resid_GE), center = T, scale = T) # note that this uses the normalized gene counts after regressing out RQN. 
summary(counts_pca)[["importance"]] %>% 
  as.data.frame() %>% 
  dplyr::select(c(PC1, PC2, PC3))  

# merge pca with metadata 
counts_pca <- as.data.frame(counts_pca$x) %>% 
  rownames_to_column(var = "LID") # value of rotated data (data * rotation matrix)
counts_pc_meta <- left_join(counts_pca, sample_info, by = "LID") # join pca values with meta data

# model PCs 
summary(lm(PC1 ~ sex + hurricane + FA_RQN + age_at_blood, counts_pc_meta)) 
summary(lm(PC2 ~ sex + hurricane + FA_RQN + age_at_blood, counts_pc_meta)) 
summary(lm(PC3 ~ sex + hurricane + FA_RQN + age_at_blood, counts_pc_meta)) 
# Hurricane is sig. for PC3, rank is sig for PC2 and 3

# Make a function that does TopGO analysis ---------------------
MW_GO_fxn <- function(my_df) {
  go_data = new('topGOdata',
                description='Simple session',
                ontology='BP',                          # biological processes 
                allGenes= my_df,                        # gene universe (all genes increasing w/ age)
                geneSelectionFun=function(x) x > 0,     # filter for sig. genes 
                nodeSize = 10,
                annotationFun = annFUN.gene2GO, 
                gene2GO = gene_ID2GO)
  
  # run fishers exact (sig/not sig. vs. GO term X/other GO term)
  go_results_fisher = runTest(go_data,
                              algorithm='weight01',
                              statistic='fisher')
  
  top_nodes <- max(go_results_fisher@geneData[4]) # I think this retrieves the number of sig. terms 
  
  # table of results 
  go_results_table <- GenTable(go_data,
                               FET.weight01 = go_results_fisher,
                               orderBy = 'FET.weight01',
                               ranksOf = 'FET.weight01',
                               topNodes = top_nodes, # Not sure why this length, could do length(age.con.go.fish@score) = 2548
                               numChar = 60)
  rownames(go_results_table)=go_results_table$GO.ID
  go_results_table$padj <- p.adjust(go_results_table$FET.weight01, method = "BH") 
  go_results_table <- go_results_table %>% 
    mutate(FET.weight01 = as.numeric(FET.weight01)) %>% 
    filter(FET.weight01 < 0.05) %>% 
    dplyr::select(c(1,2,3,4,6,7)) 
  return(go_results_table)
}

# Age: GO for age-associated genes ------------------
age_gene_univ <- emma_pre_post %>% 
  mutate(pos_age = ifelse(stdbeta_age_at_blood > 0 & qvals_age <= 0.1, 1, 0)) %>% 
  mutate(neg_age = ifelse(stdbeta_age_at_blood < 0 & qvals_age <= 0.1, 1, 0))
age_up_gene_list <- age_gene_univ$pos_age
names(age_up_gene_list) = rownames(emma_pre_post) 
s.a_age_up_GO1 <- MW_GO_fxn(age_up_gene_list) 
colnames(s.a_age_up_GO1)=c("GO.ID","Term", "annotated", "significant", "FET.weight", "enriched_pos_age_padj")

age_down_gene_list <- age_gene_univ$neg_age 
names(age_down_gene_list) = rownames(emma_pre_post) # with a lower p-val threshold there are slightly different and more sig. results 
s.a_age_down_GO1 <- MW_GO_fxn(age_down_gene_list) 
colnames(s.a_age_down_GO1)=c("GO.ID","Term", "annotated", "significant", "FET.weight", "enriched_pos_age_padj")


# Hurricane: GO for hurricane-associated genes ----------------
hurr_gene_univ <- emma_pre_post %>% 
  mutate(pos_hurricane = ifelse(stdbeta_hurricane1 > 0 & qvals_hurr <= 0.1, 1, 0)) %>% 
  mutate(neg_hurricane = ifelse(stdbeta_hurricane1 < 0 & qvals_hurr <= 0.1, 1, 0))
hurr_up_gene_list <- hurr_gene_univ$pos_hurricane
names(hurr_up_gene_list) = rownames(emma_pre_post)
hurr_up_GO1 <- MW_GO_fxn(hurr_up_gene_list) 
colnames(hurr_up_GO1)=c("GO.ID","Term", "annotated", "significant", "FET.weight", "enriched_pos_hurr_padj")

hurr_down_gene_list <- hurr_gene_univ$neg_hurricane
names(hurr_down_gene_list) = rownames(emma_pre_post)
hurr_down_GO1 <- MW_GO_fxn(hurr_down_gene_list) 
colnames(hurr_down_GO1)=c("GO.ID","Term", "annotated", "significant", "FET.weight", "enriched_pos_hurr_padj")


# Hurricane: CILP for hurricane-associated genes ---------------
cilp_hurr_results <- readRDS(file = paste0(data_dir, "CILP/cayoCILP_hurr_out.RDS"))

# clean data 
colnames(cilp_hurr_results) <- c(unlist(cilp_hurr_results[1,])[1:4], "gene1", "gene2") #set colnames (they are first line)
cilp_hurr_results <- cilp_hurr_results %>% filter(!cor_pre == "cor_pre") # filter out all colnames from concat.
cilp_hurr_results <- cilp_hurr_results %>% 
  mutate(cor_pre = as.numeric(cor_pre)) %>% 
  mutate(cor_post = as.numeric(cor_post)) %>% 
  #mutate(cor_age = as.numeric(cor_age)) %>%
  mutate(p_value = as.numeric(p_value)) %>% 
  mutate(emmreml_pval = as.numeric(emmreml_pval))

# pvals
hist(cilp_hurr_results$p_value)

# qvalue correction 
cilp_hurr_results$qval <- qvalue(p = cilp_hurr_results$emmreml_pval)$qvalues 
cilp_hurr_results <- cilp_hurr_results %>% 
  mutate(qval_low = ifelse(qval <= 0.1, "1", "0"))

### Isolate genes that sig. change corr. and find short gene name (use human analog b/c it is the same as mmul_GO but contains one more gene name)
cilp_sig <- cilp_hurr_results[cilp_hurr_results$qval < 0.1,]
cilp_sig<-left_join(cilp_sig, mmul_GO_short_unique, #mmul.hsap
                    by = c('gene1'='ensembl_gene_id')) %>% 
  dplyr::rename(hsap_homolog_gene1=external_gene_name) #hsapiens_homolog_associated_gene_name
cilp_sig$external_gene_name <- NULL
cilp_sig<-left_join(cilp_sig, mmul_GO_short_unique, 
                    by = c('gene2'='ensembl_gene_id')) %>% 
  dplyr::rename(hsap_homolog_gene2=external_gene_name)
cilp_sig[,c(3,4,8)] <- NULL #cilp_sig[,c(3,4,8,10)] (human)

mmul_GO_hurr <- mmul_GO[mmul_GO$ensembl_gene_id %in% emma_pre_post[emma_pre_post$qvals_hurr < 0.1,]$gene,]
tmp_hsps <- c(mmul_GO_hurr[str_detect(mmul_GO_hurr$name_1006, pattern = "heat shock"),]$external_gene_name, unique(mmul_GO_hurr[str_detect(mmul_GO_hurr$external_gene_name, pattern = "HSP"),]$external_gene_name), unique(mmul_GO_hurr[str_detect(mmul_GO_hurr$external_gene_name, pattern = "DNAJ"),]$external_gene_name))
binom.test(5,11,p= length(tmp_hsps)/length(unique(mmul_GO_hurr$external_gene_name)))


# Age+hurricane genes: is overlap different from background rate? ---------------------
# Fishers exact test: question = are sig. genes for aging and hurr dependent/related? 
ge_table <- rbind(c(dim(subset(emma_pre_post, qvals_age <= 0.1 & qvals_hurr <= 0.1))[1], 
                    dim(subset(emma_pre_post, qvals_age > 0.1 & qvals_hurr <= 0.1))[1]), 
                  c(dim(subset(emma_pre_post, qvals_age <= 0.1 & qvals_hurr > 0.1))[1], 
                    dim(subset(emma_pre_post, qvals_age > 0.1 & qvals_hurr > 0.1))[1]))
chisq.test(ge_table)$expected

fisher.test(ge_table)
rm(ge_table)


# Age+hurricane effects -------------------
emma_pre_post %>% 
  mutate(qvals_age_hurr = ifelse(qvals_hurr < 0.1 & qvals_age < 0.1, "low_qval_age", "high_qval_age")) %>%
  ggplot(aes(x = stdbeta_hurricane1, y = stdbeta_age_at_blood)) +
  geom_point(alpha = 0.6, size = 2.5, aes(color = qvals_age_hurr)) + #, color = wes_palette("Darjeeling2", n = 1)
  tm3_theme + 
  geom_hline(yintercept = 0, lty = 2, color = "grey70") + 
  geom_vline(xintercept = 0, lty = 2, color = "grey70") + 
  #labs(x = "Standardized effect size of hurricane", y = "Standardized effect size of age") + 
  labs(x = "", y = "") + 
  scale_x_continuous(limits = c(-11.8, 8.7)) +
  scale_y_continuous(limits = c(-11.8, 8.7)) +
  coord_equal() + 
  scale_color_manual(values = c(color_not_sig, color_both)) + 
  theme(legend.position = "none")

cor.test(x = emma_pre_post$stdbeta_hurricane1, y = emma_pre_post$stdbeta_age_at_blood)
cor.test(x = emma_pre_post$stdbeta_hurricane1, y = emma_pre_post$stdbeta_age_at_blood)$p.value

# correlation just among sig. genes
emma_pre_post_SIG <- emma_pre_post %>% 
  filter(qvals_age <= 0.1 & qvals_hurr <= 0.1)
cor.test(x = emma_pre_post_SIG$stdbeta_hurricane1, y = emma_pre_post_SIG$stdbeta_age_at_blood)

# Age+hurricane genes concordance ----------------
# are more genes concordant than expected? Yes 
emma_pre_post_SIG <- emma_pre_post %>% 
  filter(qvals_age <= 0.1 & qvals_hurr <= 0.1)

ge_table <- rbind(c(dim(subset(emma_pre_post_SIG, 
                               stdbeta_age_at_blood > 0 & stdbeta_hurricane1 > 0))[1], 
                    dim(subset(emma_pre_post_SIG, 
                               stdbeta_age_at_blood > 0 & stdbeta_hurricane1 < 0))[1]), 
                  c(dim(subset(emma_pre_post_SIG, 
                               stdbeta_age_at_blood < 0 & stdbeta_hurricane1 < 0))[1], 
                    dim(subset(emma_pre_post_SIG, 
                               stdbeta_age_at_blood < 0 & stdbeta_hurricane1 > 1))[1]))
chisq.test(ge_table)$expected

fisher.test(ge_table)
rm(ge_table, emma_pre_post_SIG)


# concordance plots 
emma_pre_post$concordance = with(emma_pre_post, stdbeta_age_at_blood*stdbeta_hurricane1 > 0) 
# Calculate concordance statistics, with 1000 bootstrap replicates to estimate confidence intervals
emma_pre_post_concordance = rbind(
  do.call(rbind,lapply(seq(0.01,1,0.01),function(x) {
    pass.filter = subset(emma_pre_post, qvals_age <= x & qvals_hurr <= x) 
    out = replicate(1000,mean(sample(pass.filter$concordance,replace=TRUE)))
    data.frame(
      qval=x,                                          # pval threshold (change to FDR**)
      concordant=with(pass.filter,sum(concordance)),   # Number of concordant genes/sites
      concordance=with(pass.filter,mean(concordance)), # Fraction of genes that are concordant
      stdev=with(pass.filter,sd(concordance)),         # Standard deviation of concordance
      n=nrow(pass.filter),                             # Number of sites at this FDR threshold
      error.min=quantile(out,0.025),                   # Minimum of 95% confidence interval
      error.max=quantile(out,0.975),                   # Maximum of 95% confidence interval
      stringsAsFactors=FALSE
    )
  }))
)

concord_tmp<-emma_pre_post_concordance %>% 
  ggplot(aes(y = concordance, x = qval)) + 
  geom_point() + 
  tm3_theme + 
  labs(x = "FDR threshold", y = "Fraction concordant") + 
  scale_y_continuous(limits = c(0.5, 1)) + 
  geom_errorbar(aes(ymin=error.min,ymax=error.max))

