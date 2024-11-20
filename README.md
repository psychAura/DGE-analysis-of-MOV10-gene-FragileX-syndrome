# Differential Gene Expression analysis of MOV 10 gene in Fragile X Syndrome

### workshop from Havard Chan Bioinformatics core.

An end-to-end gene-level RNA-seq differential expression workflow using various R packages. Workflow:

- Read data obtained from Salmon, 

-convert pseudocounts to counts, 

-perform exploratory data analysis for quality assessment and to explore the relationship between samples,

- perform differential expression analysis, 

-and visually explore the results prior to performing downstream functional analysis.

## Dataset

 The data is a publicly available RNA-Seq dataset that is part of a larger study described in Kenny PJ et al, Cell Rep 2014.

The RNA-Seq was performed on HEK293F cells that were either transfected with a MOV10 transgene, or siRNA to knock down Mov10 expression, or non-specific (irrelevant) siRNA. This resulted in 3 conditions Mov10 oe (over expression), Mov10 kd (knock down) and Irrelevant kd, respectively.

Using these data, we will evaluate transcriptional patterns associated with perturbation of MOV10 expression. Please note that the irrelevant siRNA will be treated as our control condition.

## Dependencies
## Setup
### Bioconductor and CRAN libraries used
-library(DESeq2)

-library(tidyverse)

-library(RColorBrewer)

-library(pheatmap)

-library(DEGreport)

-library(tximport)

-library(ggplot2)

-library(ggrepel)

## Primary folder

-de_scrippt.R