## Gene-level differential expression analysis using DESeq2
#load libraries
library(DESeq2)
library(tidyverse)
library(RColorBrewer)
library(pheatmap)
library(ggrepel)
library(cowplot)
library(clusterProfiler)
library(DEGreport)
library(org.Hs.eg.db)
library(DOSE)
library(pathview)
library(tximport)
library(AnnotationHub)
library(ensembldb)
library(ggplot2)

sessionInfo()

## List all directories containing data  
samples <- list.files(path = "./data", full.names = T, pattern="salmon$")

## Obtain a vector of all filenames including the path
files <- file.path(samples, "quant.sf")

## Since all quant files have the same name it is useful to have names for each element
names(files) <- str_replace(samples, "./data/", "") %>% 
  str_replace(".salmon", "")

# Load the annotation table for GrCh38
tx2gene <- read.delim("tx2gene_grch38_ens94.txt")

# Take a look at it 
tx2gene %>% View()

?tximport   # let's take a look at the arguments for the tximport function

# Run tximport
txi <- tximport(files, type="salmon", tx2gene=tx2gene[,c("tx_id", "ensgene")], countsFromAbundance="lengthScaledTPM")

attributes(txi)

# Look at the counts
txi$counts %>% View()

# Write the counts to an object
data <- txi$counts %>% 
  round() %>% 
  data.frame() 

## Create a sampletable/metadata
sampletype <- factor(c(rep("control",3), rep("MOV10_knockdown", 2), rep("MOV10_overexpression", 3)))
meta <- data.frame(sampletype, row.names = colnames(txi$counts))

ggplot(data) +
  geom_histogram(aes(x = Mov10_oe_1), stat = "bin", bins = 200) +
  xlab("Raw expression counts") +
  ylab("Number of genes")

mean_counts <- apply(data[,6:8], 1, mean)        #The second argument '1' of 'apply' function indicates the function being applied to rows. Use '2' if applied to columns 
variance_counts <- apply(data[,6:8], 1, var)
df <- data.frame(mean_counts, variance_counts)

ggplot(df) +
  geom_point(aes(x=mean_counts, y=variance_counts)) + 
  scale_y_log10(limits = c(1,1e9)) +
  scale_x_log10(limits = c(1,1e9)) +
  geom_abline(intercept = 0, slope = 1, color="red")

### Check that sample names match in both files
all(colnames(txi$counts) %in% rownames(meta))
all(colnames(txi$counts) == rownames(meta))

## Create DESeq2Dataset object
dds <- DESeqDataSetFromTximport(txi, colData = meta, design = ~ sampletype)

##plotting PCA without DesEq
dds <- estimateSizeFactors(dds)
sizeFactors(dds)

normalized_counts <- counts(dds, normalized=TRUE)
write.table(normalized_counts, file="data/normalized_counts.txt", sep="\t", quote=F, col.names=NA)

### Transform counts for data visualization
rld <- rlog(dds, blind=TRUE)

### Plot PCA 
plotPCA(rld, intgroup="sampletype")

###hierarchical clustering for mov10
### Extract the rlog matrix from the object
rld_mat <- assay(rld)    

## "assay()" is part of the "SummarizedExperiment" package which is a DESeq2 dependency and is loaded with the DESeq2 library
### Compute pairwise correlation values
rld_cor <- cor(rld_mat)    ## cor() is a base R function

## check the output of cor(), make note of the row names and column names
head(rld_cor)

head(meta)
### Load pheatmap package
library(pheatmap)

### Plot heatmap using the correlation matrix and the metadata object
pheatmap(rld_cor, annotation = meta)

dds <- DESeqDataSetFromTximport(txi, colData = meta, design = ~ sampletype)
dds <- DESeq(dds)

## Total number of raw counts per sample
colSums(counts(dds))

## Total number of normalized counts per sample
colSums(counts(dds, normalized=T))

## Plot dispersion estimates
plotDispEsts(dds)

## Define contrasts for MOV10 overexpression
contrast_oe <- c("sampletype", "MOV10_overexpression", "control")

## Extract results for MOV10 overexpression vs control
res_tableOE <- results(dds, contrast=contrast_oe, alpha = 0.05)

res_tableOE %>% 
  data.frame() %>% 
  View()

## Save the unshrunken results to compare
res_tableOE_unshrunken <- res_tableOE

# Apply fold change shrinkage
res_tableOE <- lfcShrink(dds, coef="sampletype_MOV10_overexpression_vs_control", type="apeglm")

# MA plot using unshrunken fold changes
plotMA(res_tableOE_unshrunken, ylim=c(-2,2))
 
# MA plot using shrunken fold changes
plotMA(res_tableOE, ylim=c(-2,2))

## Define contrasts for knockdown samples
contrast_kd <- c("sampletype", "MOV10_knockdown", "control")

## Extract results for MOV10 knockdown vs control
res_tableKD <- results(dds, contrast=contrast_kd, alpha = 0.05)

res_tableKD %>% 
  data.frame() %>% 
  View()

## Save the unshrunken results to compare
res_tableKD_unshrunken <- res_tableKD

# Apply fold change shrinkage
res_tableKD <- lfcShrink(dds, coef="sampletype_MOV10_knockdown_vs_control", type="apeglm")

# MA plot using unshrunken fold changes
plotMA(res_tableKD_unshrunken, ylim=c(-2,2))

# MA plot using shrunken fold changes
plotMA(res_tableKD, ylim=c(-2,2))

## Summarize results for MOV10 OE
summary(res_tableOE, alpha = 0.05)

##extracting significant diffrentiall exressed genes
####set thresholds
padj.cutoff <- .05

# Create a tibble of results
res_tableOE_tb <- res_tableOE %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

# Subset the tibble to keep only significant genes
sigOE <- res_tableOE_tb %>%
  dplyr::filter(padj < padj.cutoff)

# Take a quick look at this tibble
sigOE

# Create a tibble of results
res_tableKD_tb <- res_tableKD %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

# Subset the tibble to keep only significant genes
sigKD <- res_tableKD_tb %>%
  dplyr::filter(padj < padj.cutoff)
