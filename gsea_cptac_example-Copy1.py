#!/usr/bin/env python
# coding: utf-8

# Gene Set Enrichment Analysis
# - GSEA on protein data obtained from standard and immunoproteasome dataframes. The aim is to identify the upregulated genes in standard and immunoproteasome samples and determine their enriched biological pathways. We will use the gseapy library for GSEA, and the usual data science libraries for data manipulation and visualization.

# Step 1: Importing packages and setting up your notebook.¶
# - The Python environment needs several libraries for data manipulation and visualization.

pip install gseapy


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import gseapy as gp
from gseapy.plot import barplot, dotplot


import os
#path="/home/eewon.tai01/Downloads/GSEAPY"
path = "gsea_standard_immunoproteasome"
os.chdir(path)


### Data format conversion - into a form gesapy would accept
# convert into pandas dataframe
# changed these paths
tmt_labels = pd.read_csv("TMT_labels.csv", dtype=object)
exp_data = pd.read_csv("cleaned_and_normalised.csv", dtype=object)

# for simplicity, take only rows that have all values
exp_data.dropna(inplace=True)


# Step 2: Data Segregation
# - To prepare the data for analysis, we split the proteomics data into two separate dataframes, one for SP samples and another for IP samples.

## Automate process of searching for column names that start with SP or IP
## Also automate process of searching for column indices that specify the above
columns = list(exp_data.columns)
class_vector_sp = list()
class_vector_ip = list()
indices_sp = list()
indices_ip = list()
i=0
for col in columns:
    col = col.split("_")[0]
    if col == "SP":
        class_vector_sp.append(col)
        indices_sp.append(i)
    elif col == "IP":
        class_vector_ip.append(col)
        indices_ip.append(i)
    i += 1
print("class_vector_sp", class_vector_sp)
print("\n\n\n\n\nclass_vector_ip", class_vector_ip)
print("\n\n\n\n\nindices_sp", indices_sp)
print("\n\n\n\n\nindices_ip", indices_ip)

data_genenames = exp_data.iloc[:, 0]
data_sp = exp_data.iloc[:, indices_sp]
data_ip = exp_data.iloc[:, indices_ip]
# The problem was my dataframe was an "object" instead of numbers so I was getting an empty dataframe for the computation. 
#data2 = data2.astype(float)
#data = pd.concat([data1, data2], axis=1)
#print(data)

print("\n\n\n\n\ndata_genenames", data_genenames)
print("\n\n\n\n\ndata_sp", data_sp)
print("\n\n\n\n\ndata_ip", data_ip)

standard = pd.concat([data_genenames, data_sp], axis=1)
immuno = pd.concat([data_genenames, data_ip], axis=1)

print("\n\n\n\n\nstandard", standard)
print("\n\n\n\n\nimmuno", immuno)


# in cptac example (use case 5) - corresponds to:
# standard - normal
# immuno - tumor


# Step 3: Performing Welch's t-test
# - To identify upregulated genes, we use Welch's t-test, a version of the two-sample t-test that accommodates different variances between two groups. In our case, we compare the gene abundance of SP and IP groups for each gene. Genes with a significant p-value after correction for multiple testing are considered upregulated, and their expression patterns (i.e., higher in SP or IP samples) are noted.

# Create array variables to hold the significant genes for each partition
immuno_genes = []
standard_genes = []

# Grab the genes of interest, ignoring the MSI column in the dataframe
genes = data_genenames

# Correct alpha level for multiple testing by dividing the standard .05 by the number of genes to be analyzed
threshold = .05 / len(genes)

# Perform Welch's t-test (different variances) on each gene between the two groups
for i in range(len(data_genenames)):
#    gene = genes[i]
    immuno_gene_abundance = pd.to_numeric(data_ip.iloc[i, ], downcast='float')
    standard_gene_abundance = pd.to_numeric(data_sp.iloc[i, ], downcast='float')

    ttest = stats.ttest_ind(immuno_gene_abundance, standard_gene_abundance, equal_var=False)
    pvalue = ttest.pvalue
    # If the P-value is significant, determine which partition is more highly expressed
    if pvalue < threshold:
        if immuno_gene_abundance.mean() > standard_gene_abundance.mean():
            immuno_genes.append(genes[i])
        elif standard_gene_abundance.mean() > immuno_gene_abundance.mean():
            standard_genes.append(genes[i])

print("Proteomics immuno Genes:", len(immuno_genes))
print("Proteomics standard Genes:", len(standard_genes))

print("\n\n\n\n\nProteomics immuno Genes:", immuno_genes)
print("\n\n\n\n\nProteomics standard Genes:", standard_genes)

print("\n\n\n\n\ntumor_gene_abundance:", immuno_gene_abundance)
print("\n\n\n\n\nnormal_gene_abundance:", standard_gene_abundance)


# Suppose we observe two independent samples, e.g. flower petal lengths, and we are considering whether the two samples were drawn from the same population (e.g. the same species of flower or two species with similar petal characteristics) or two different populations.
# 
# The t-test quantifies the difference between the arithmetic means of the two samples. The p-value quantifies the probability of observing as or more extreme values assuming the null hypothesis, that the samples are drawn from populations with the same population means, is true. A p-value larger than a chosen threshold (e.g. 5% or 1%) indicates that our observation is not so unlikely to have occurred by chance. Therefore, we do not reject the null hypothesis of equal population means. If the p-value is smaller than our threshold, then we have evidence against the null hypothesis of equal population means.

print("ttest.statistic", ttest.statistic)  # t-statistic
print("ttest.df", ttest.df)  # The number of degrees of freedom used in calculation of the t-statistic
print(ttest.confidence_interval(confidence_level=0.95))
# Computes a confidence interval around the difference in population means for the given confidence level.


# Step 4: Gene set enrichment analysis (GSEA)
# - With our list of upregulated genes, we perform a GSEA using the enrichr() function from the gseapy library. This analysis identifies biological pathways that are overrepresented in our list of genes, providing insight into potential molecular mechanisms at play.

immuno_enr = gp.enrichr(gene_list=immuno_genes, gene_sets='ENCODE_ChIP-seq.gmt', outdir='gsea_standard_immunoproteasome/enrichr_kegg_immuno')
standard_enr = gp.enrichr(gene_list=standard_genes, gene_sets='ENCODE_ChIP-seq.gmt', outdir='gsea_standard_immunoproteasome/enrichr_kegg_standard')

# View the data as a table
print(immuno_enr.res2d.head())
print(standard_enr.res2d.head())


# Step 5: Visualizing Enrichment Results
# - We create barplots to visualize the GSEA results, giving us a better understanding of the pathways significantly associated with the upregulated genes in both the SP and IP groups.

barplot(immuno_enr.res2d, title="Proteomics Immuno KEGG_2016")
plt.show()

barplot(standard_enr.res2d, title="Standard KEGG_2016")
plt.show()


# Discussion
# - https://maayanlab.cloud/Enrichr/enrich
# - https://maayanlab.cloud/chea3/
# - https://www.gsea-msigdb.org/gsea/msigdb/mouse/genesets.jsp
# - https://en.wikipedia.org/wiki/Welch%27s_t-test
# - https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.ttest_ind.html



