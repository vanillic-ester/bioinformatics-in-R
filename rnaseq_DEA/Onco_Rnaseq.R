# Set working directory
setwd("~/Desktop/rnaseq_DEA") # this should be set to where the input files are located

# necessary libraries
library(dplyr)
library(tidyverse)
library(readr)
library(GEOquery)
library(DESeq2)
library(EnhancedVolcano)

# -----------------------------
# Load and Preprocess the Data
# -----------------------------

# read the data
data <- read.csv(file = 'GSE287307_raw_counts_all_samples.txt', sep = "\t")

# get metadata
gse <- getGEO(GEO='GSE287307',GSEMatrix = TRUE)
metadata <- pData(phenoData(gse[[1]]))

metadata_subset <- select(metadata,c(1,43)) # contains sample and genotype

col_data <- metadata_subset
colnames(col_data)[2] <- "genotype"

# extract read counts from data
count_data <- select(data,c(1:7))

# format counts for differential expression
rownames(count_data) <- count_data[,1]
count_data <- count_data[,-1]

# change row names for col_data to match our columns of counts
rownames(col_data) <- colnames(count_data)

# make sure that rows of col_data and columns of counts match
all(rownames(col_data) == colnames(count_data))

# -----------------------------
# Perform Differential Gene expression Analysis and Modeling 
# -----------------------------

# make DEseq2 object
dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData= col_data,
                              design = ~ genotype)

# prefiltering steps

filtered_dds <- dds[rowSums(counts(dds)) > 10, ]

# set factor level
filtered_dds$genotype <- relevel(filtered_dds$genotype, ref='WT')

# run DESeq (differential gene expression)
dds_out <- DESeq(filtered_dds)

# extract results and order them to identify differentially expressed genes
res <- results(dds_out)
res <- res[order(res$padj),]

# -----------------------------
# Visualize differential gene expression
# -----------------------------

# Create a mapping of Ensembl IDs to Gene Names
id_to_names <- setNames(data$gene_name, data$gene_id)

### Data Visualization
# necessary libraries for visualization
library(genefilter)
library(pheatmap)
library(org.Mm.eg.db)

# transform count data for visualization
rld <- rlogTransformation(dds_out, blind=F)

# extract top genes for visualization
top_genes <- head(order(rowVars(assay(rld)), decreasing = TRUE), 20)
mat <- assay(rld)[top_genes, ]

# convert rows to gene names to make heatmap readable
rownames(mat) <- id_to_names[rownames(mat)]

# use heatmap to visualize gene expression
pheatmap(mat, scale='row')
# use volcano plot to visualize significant differential expressed genes
EnhancedVolcano(res,lab = rownames(res), x='log2FoldChange', y='pvalue')

# -----------------------------
# Functional Annotation and Enrichment Analysis
# -----------------------------

# Load necessary libraries for enrichment analysis
library(clusterProfiler)
library(DOSE)

# filter for upregulated and downregulated knockout genes 
up <- subset(res,res$padj < 0.05 & res$log2FoldChange > 0)
down <- subset(res,res$padj < 0.05 & res$log2FoldChange < 0)

# output files of the gene ids which can be used for GO analysis in DAVID
write(data$gene_id, "background.txt")
write(rownames(down),"down.txt")
write(rownames(up),"up.txt")

# Perform Gene Ontology (GO) enrichment analysis
up_enriched <- enrichGO(rownames(up), OrgDb = org.Mm.eg.db, ont = 'ALL', pAdjustMethod = 'BH', pvalueCutoff = 0.05,
                        qvalueCutoff = 0.2, keyType = 'ENSEMBL')

down_enriched <- enrichGO(rownames(down), OrgDb = org.Mm.eg.db, ont = 'ALL', pAdjustMethod = 'BH', pvalueCutoff = 0.05,
                        qvalueCutoff = 0.2, keyType = 'ENSEMBL')

# view enriched pathways from upregulated and downregulated genes from knockouts
barplot(down_enriched, showCategory=8, drop=T)
barplot(up_enriched, showCategory=8, drop=T)
