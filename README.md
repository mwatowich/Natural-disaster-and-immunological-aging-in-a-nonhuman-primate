# Natural Disaster and Immunological Aging in a Nonhuman Primate

This repository contains code written for an investigation of the effects of aging and Hurricane Maria on rhesus macaque gene expression from the peripheral blood immune system.

> Watowich MM, Chiou KL, Montague MJ, Simons ND, Horvath JE, Ruiz-Lambides AV, Martínez MI, Higham JP, Brent LJN, Platt ML, Snyder-Mackler N. Natural disaster and immunological aging in a nonhuman primate. In Press. 2022.


## Pipeline

1. First download and index the reference genome. Map reads and quantify expression.

	Required software: [Kallisto](https://pachterlab.github.io/kallisto) (v0.43.1)

```
scripts/kallisto_mapping.sh
```

2. Import expression matrix and filter and normalize counts.

	Required software: [R](https://cran.r-project.org)

```
scripts/normalize_gene_counts.R
```

3. Perform linear modeling at the gene-level and perform enrichment tests.

	Required software: [R](https://cran.r-project.org)

```
scripts/rnaseq_linear_modeling.R
scripts/enrichment_tests.R
```

4. Estimate transcriptomic age from human data.

	[Please see] (https://github.com/smacklab/cayo_santiago_immune_aging_review/blob/master/scripts/rna_human_comparison.R)


5. Finally, use single cell data generated from PBMCs in rhesus macaques to identify top marker genes of immune cells.
