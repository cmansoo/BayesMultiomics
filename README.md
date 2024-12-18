# BayesMultiomics
 BayesMultiomics is designed for statistical modeling and inference with a two stage Bayesian shrinkage model specially designed for high dimensional multi-platform genomics and health care data. 
  Our R package BayesMultiomics is based on Bayesian shrinkage based data integration model developed by us earlier (Xue et al., 2024). 
  This model can integrate and synthesize data from four different data platforms as a) gene expression, b) DNA Methylation, c) gene functional classification obtained from Database for Annotation, Visualiza-tion, and Integrated Discovery (DAVID) (Sherman et al., 2022; Huang et al., 2009), and d) clinical features and clinical outcomes. 
  Our BayesMultiomics package offers versatility and flexibility by incorporating information across all 4 mentioned platforms to learn the underlying biological mechanisms among different high-throughput platforms. 
  In addition, our package is scalable and can perform high-dimensional variable selection for discovering relationship across platforms and identifying clinically relevant biomarkers.

## Installation
```
remotes::install_github("cmansoo/BayesMultiomics", build_vignettes=TRUE)
```

## Vignettes
```
library(BayesMultiomics)
vignette("BayesMultiomics")
vignette("DAVIDGeneGrouping")
```
