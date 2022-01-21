### Enrichment analyses

# Load metadata
metadata<-as.data.frame(read_xlsx('~/Downloads/PNAS_202121663_s2.xlsx', skip = 1, col_names = T))

# Load model results
rna_mod_out<-as.data.frame(read_xlsx('~/Downloads/PNAS_202121663_s3.xlsx', skip = 1, col_names = T))


## Fishers Exact Test of overlap between age and hurricane exposure effects
ge_table<-rbind(c(dim(subset(rna_mod_out, `FDR aging` <= 0.1 & `FDR hurricane exposure` <= 0.1))[1], 
                  dim(subset(rna_mod_out, `FDR aging` > 0.1 & `FDR hurricane exposure` <= 0.1))[1]), 
                c(dim(subset(rna_mod_out, `FDR aging` <= 0.1 & `FDR hurricane exposure` > 0.1))[1], 
                  dim(subset(rna_mod_out, `FDR aging` > 0.1 & `FDR hurricane exposure` > 0.1))[1]))
chisq.test(ge_table)$expected
fisher.test(ge_table)

## Correlation of standardized betas for aging and hurricane effects 
cor.test(rna_mod_out$`standardized beta aging`, rna_mod_out$`standardized beta hurricane exposure`)
cor.test(rna_mod_out[rna_mod_out$`FDR aging`<0.1 & 
                       rna_mod_out$`FDR hurricane exposure`<0.1,]$`standardized beta aging`, 
         rna_mod_out[rna_mod_out$`FDR aging`<0.1 & 
                       rna_mod_out$`FDR hurricane exposure`<0.1,]$`standardized beta hurricane exposure`)

## Binomial test of direction of age and hurricane effects at the gene-level
binom.test(dim(subset(rna_mod_out, `FDR aging` <= 0.1 & `FDR hurricane exposure` <= 0.1 & 
                        (`beta aging` > 0 & `beta hurricane exposure` > 0 | 
                           `beta aging` < 0 & `beta hurricane exposure` < 0)))[1], 
           dim(subset(rna_mod_out, `FDR aging` <= 0.1 & `FDR hurricane exposure` <= 0.1))[1], p = 0.5)


### Gene Ontology functional enrichment tests --------------------
## Example: genes that significantly increase in expression with aging
# Make list of genes that increase with aging and are significantly associated with aging
go_age_genes <- rna_mod_out %>% 
  mutate(pos_age = ifelse(`standardized beta aging` > 0 & `FDR aging` <= 0.1, 1, 0))
go_age_genes_v <- go_age_genes$pos_age
names(go_age_genes_v) = rna_mod_out$`ensembl ID`

# Create biomaRt object
library(biomaRt)
mmul_mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                              dataset = "mmulatta_gene_ensembl", 
                              host = "www.ensembl.org")
# Create biomaRt object with each genes' associated GO terms
mmul_GO <- biomaRt::getBM(attributes = c('ensembl_gene_id','go_id','name_1006','external_gene_name'),
                          filters = 'ensembl_gene_id',
                          values = row.names(norm_gene_counts), 
                          mart = mmul_mart,
                          uniqueRows = TRUE)

# Make list of genes and associated GO IDs
gene_ID2GO <- lapply(unique(mmul_GO$ensembl_gene_id),
                     function(x){sort(mmul_GO[mmul_GO$ensembl_gene_id == x,
                                              'go_id'])})
names(gene_ID2GO) <- unique(mmul_GO$ensembl_gene_id)

# Gene Ontology enrichment analysis
library(topGO)
go_data = new(Class = 'topGOdata',
              description='Simple session',
              ontology='BP',
              allGenes= go_age_genes_v,
              geneSelectionFun=function(x) x > 0,
              nodeSize = 10,
              annotationFun = annFUN.gene2GO, 
              gene2GO = gene_ID2GO)
go_results_fisher<-runTest(go_data,
                           algorithm='weight01',
                           statistic='fisher')
top_nodes <- max(go_results_fisher@geneData[4])
go_results_table <- GenTable(go_data,
                             FET.weight01 = go_results_fisher,
                             orderBy = 'FET.weight01',
                             ranksOf = 'FET.weight01',
                             topNodes = top_nodes,
                             numChar = 60)
go_results_table$padj <- p.adjust(go_results_table$FET.weight01, method = "BH")



## Transcription factor binding site enrichment tests ---------------------
# Example: genes that significantly increase with aging
tfbs_gene_list<-data.frame(gene = rna_mod_out[rna_mod_out$`FDR aging`<0.1 & 
                                                rna_mod_out$`standardized beta aging`>0,]$`ensembl ID`)
tfbs_bkgd<-data.frame(subset(rna_mod_out, !(`ensembl ID` %in% homer_gene_list$gene), `ensembl ID`))
# Run HOMER on the command line with standard inputs

