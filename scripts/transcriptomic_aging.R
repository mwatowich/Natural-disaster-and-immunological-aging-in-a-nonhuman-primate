
peters_sample_info <- sample_info[sample_info$age_at_blood > 1,] 
peters_resid_GE <- resid_GE[,colnames(resid_GE) %in% 
                              peters_sample_info$LID]

# scale the expression data 
e.peters = t(apply(peters_resid_GE,1,scale)) # scale the normalized counts w/in gene, all individuals 
dimnames(e.peters) = list(rownames(peters_resid_GE),colnames(peters_resid_GE)) 

# read in predictors from Peters et al. and filter to include genes in the rhesus dataset 
peters.predict = read.delim(paste0(data_dir,'Peters_transcriptional_aging/peters_predictors.txt'))
peters.predict = merge(peters.predict,mmul.hsap,by.x='gene_id',by.y='hsapiens_homolog_associated_gene_name')
peters.predict = subset(peters.predict,ensembl_gene_id %in% rownames(e.peters))

# get rid of genes that have multiple macaque ensembl IDs (first, keep all the genes with 1 copy, and then select just a single row of the genes that are represented multiple times). 
peters.predict = rbind(
  subset(peters.predict,!ensembl_gene_id %in% names(which(table(peters.predict$ensembl_gene_id)>1))),
  do.call(rbind,lapply(split(subset(peters.predict,ensembl_gene_id %in% names(which(table(peters.predict$ensembl_gene_id)>1))),subset(peters.predict,ensembl_gene_id %in% names(which(table(peters.predict$ensembl_gene_id)>=2)))$ensembl_gene_id),function(x) x[1,]))
) ## MMW changed ==2 to >=2 b/c there was a warning b/c one ensembl_gene_id had 3 assoc. genes 
rownames(peters.predict) = peters.predict$ensembl_gene_id

# get only the genes which were age-related with pre-/post- model 
age_assoc_genes <- rownames(subset(emma_pre_post, qvals_age < 0.1)) 
peters.predict = peters.predict[rownames(peters.predict) %in% age_assoc_genes,] # age_assoc_genes is selected by a qval < 0.1. Can use either emma_pre_post or age_hurr_coeffs 
e.peters = e.peters[rownames(e.peters) %in% rownames(peters.predict), colnames(e.peters) %in% peters_sample_info$LID] # genes in Peters and our study, and only individuals from the sub-sample - try all animals vs. subsample (w/ less 5yr) vs. adults 

# Equalize gene orders 
e.peters = e.peters[sort(rownames(e.peters)),]
peters.predict = peters.predict[rownames(e.peters),] 

# Save predictions, bind age.zs and real age to plot
peters.predictions = subset(peters_sample_info,select=c('animal_id','LID','RID','sex','age_at_blood', 'hurricane', 'FA_RQN')) 

# Perform prediction across a range of FDR thresholds (K Chiou code)
# age_coeffs is a data frame with model outputs for macaque datasets. So the effect sizes (Z and betas), p values, q values (FDRs), etc were all in there for each orthologous gene. 
peters.predictions.fdr = lapply(seq(0.01,0.1,0.01),function(x) { #can go to higher qvals but we have already filtered for age-associated genes that overlap with Peters (above) so it doesnt make sense to do so
  e.peters.fdr = e.peters[rownames(e.peters) %in% rownames(subset(emma_pre_post,qvals_age <= x)),] # finding the age associated genes. either use age_hurr_coeffs of emma_pre_post 
  peters.predict.fdr = peters.predict[rownames(e.peters.fdr),]
  age.z.fdr = unlist(lapply(colnames(e.peters.fdr),function(x) sum(e.peters.fdr[,x] * peters.predict.fdr$z)))
  names(age.z.fdr) = colnames(e.peters.fdr)
  age.zs.fdr = mean(peters_sample_info$age_at_blood) + (age.z.fdr - mean(age.z.fdr)) * sd(peters_sample_info$age_at_blood) / sd(age.z.fdr) 
  out = peters.predictions
  out$qval = x
  out$zs = age.zs.fdr[out$LID]
  out
})
# mean chron age + deviation of individual's predicted age from mean predicted age * sd(chron age)/sd(bio age)
# Calculate fit across FDR ranges
peters.predictions.fit = do.call(rbind,lapply(peters.predictions.fdr,function(x) {
  m = lm(zs~age_at_blood,data=x)
  data.frame(qval=unique(x$qval),n=sum(with(emma_pre_post,qvals_age < unique(x$qval))),b=coef(m)[2],r2=summary(m)$r.squared) # either use emma_pre_post or age_hurr_coeffs
}))

# Plot Rsquared and b of fit 
ggplot(peters.predictions.fit,aes(qval,r2)) + geom_point() + theme_classic() + scale_y_continuous(limits=c(0,1)) 
ggplot(peters.predictions.fit,aes(qval,b)) + geom_point() + theme_classic() + scale_y_continuous(limits=c(0,1))

# Bind results from all FDRs, select FDR threshold we want 
peters.predictions.all = do.call(rbind,peters.predictions.fdr) %>% 
  mutate(age_diff = (zs-age_at_blood))
peters.predictions.using <- peters.predictions.all[peters.predictions.all$qval == "0.01",] # stringent age-assoc genes 

# Find correlation between age predictions and chronological age - for pre-hurricane samples (to see how good the predictions are)
##cor.test(peters.predictions.using$age_at_blood, peters.predictions.using$zs)#$p.value

# Plot chronological age vs. predicted age
# age calculations change with different fdr cut-offs 
peters.predictions.using %>% 
  ggplot(aes(x = age_at_blood, y = zs, color = hurricane)) + 
  geom_point(alpha = 0.5, size = 2) +
  geom_smooth(method=lm, show.legend = F) + 
  geom_smooth(method=lm, se = F) + 
  tm3_theme + 
  scale_color_manual(values = wes_palette("Darjeeling1", n = 5)[c(4,5)], 
                     labels = c("Pre-\nhurricane", "Post-\nhurricane"), name = "") + 
  labs(x = "Chronological age", y = "Biological age") + 
  theme(legend.key.height=unit(4,"line"))

# model pred_age ~ chron+hurr, subtract out chron_age partial residuals to get effect of hurr on pred_age, plot age diff due to hurr. 
mod_pred_age <- as.data.frame(summary(lm(peters.predictions.using$zs ~ peters.predictions.using$age_at_blood + peters.predictions.using$hurricane))$coefficients) # sex wasnt sig. either. 
peters.predictions.using$resid_zs <- peters.predictions.using$zs - (mod_pred_age[2,1]*peters.predictions.using$age_at_blood) - mod_pred_age[1,1] # resid_age (due to hurr) = bio age - beta(age)*chron age - intercept 
mean(peters.predictions.using[peters.predictions.using$hurricane == 1, "resid_zs"]) - 
  mean(peters.predictions.using[peters.predictions.using$hurricane == 0, "resid_zs"])

# plot effect of hurr on partial residuals of predicted age 
pre_predict_age_median <- median(peters.predictions.using[peters.predictions.using$hurricane == "0", "resid_zs"])
post_predict_age_median <- median(peters.predictions.using[peters.predictions.using$hurricane == "1", "resid_zs"])
peters_tmp<-peters.predictions.using %>% 
  arrange(hurricane) %>%
  ggplot(aes(x = resid_zs, y = hurricane, color = hurricane, fill = hurricane)) + 
  geom_density_ridges(alpha = 0.8, scale = 1.5, bandwidth = 2, show.legend = F) + #scale=2
  tm3_theme + 
  geom_segment(aes(x=pre_predict_age_median, 
                   xend = pre_predict_age_median, y=0.55, yend=4), #y=0.55, yend=4
               color = color_pre, lty = 2, size = 2) + #wes_palette("Darjeeling1", n = 5)[c(4)]
  geom_segment(aes(x=post_predict_age_median, 
                   xend = post_predict_age_median, y=0.55, yend=4), #y=0.55, yend=4
               color = color_post, lty = 2, size = 2) + #wes_palette("Darjeeling1", n = 5)[c(5)]
  scale_color_manual(values = c(color_pre, color_post)) + #wes_palette("Darjeeling1", n = 5)[c(4,5)]
  scale_fill_manual(values = c(color_pre, color_post)) + #wes_palette("Darjeeling1", n = 5)[c(4,5)]
  labs(x = "Difference between biological and chronological age", y = "Relation to Hurricane Maria") + 
  scale_y_discrete(labels = c("pre-", "post-")) +
  scale_x_continuous(labels = c("0"), breaks = 0) +
  theme(axis.title.y = element_text(angle = 90, hjust = 0.4, size = 28),
        axis.ticks.y = element_blank(), 
        axis.title.x = element_text(size = 24),
        axis.text.y = element_text(size = 24))
