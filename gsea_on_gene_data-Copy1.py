#!/usr/bin/env python
# coding: utf-8

pip install gseapy


#pip install gseapy
# install gseapy package


import gseapy as gp

names = gp.get_library_name()
print(names)
# print library names already in gseapy package
# we are interested in GO (Gene Ontology) gene sets:
# 'GO_Biological_Process_2023', 'GO_Cellular_Component_2023', 'GO_Molecular_Function_2023'
# we focus on Biological Process, because we are interested in what processes change when Standard proteasome becomes immunoproteasome


# import packages
import pandas as pd
import numpy as np
import gseapy as gp
import matplotlib.pyplot as plt
from gseapy.plot import barplot, dotplot

### Data format conversion - into a form gesapy would accept
# convert into pandas dataframe
tmt_labels = pd.read_csv("TMT_labels.csv", dtype=object)
exp_data = pd.read_csv("cleaned_and_normalised.csv", dtype=object)

# for simplicity, take only rows that have all values
exp_data.dropna(inplace=True)

## Automate process of searching for column names that start with SP or IP
## Also automate process of searching for column indices that specify the above
# SP=Standard Proteasome, IP=Immunoproteasome
columns = list(exp_data.columns)
class_vector = list()  # vector that contains the class values for all the columns that are SP or IP - either SP or IP
indices = list()   # indices of columns of interest (columns that are either SP or IP)
i=0
for col in columns:
    col = col.split("_")[0]
    if col == "SP" or col == "IP":  # if column starts with "SP" or "IP", append it to class_vector and indices
        class_vector.append(col)
        indices.append(i)
    i += 1
print("class_vector: \n", class_vector)
print("indices: \n", indices)

data1 = exp_data.iloc[:, 0] # all gene names
print("data1 (gene names): \n", data1)

data2 = exp_data.iloc[:, indices]  # dataframe of columns that we are interested in

# The dataframe is an "object" instead of numbers so I was getting an empty dataframe
# type conversion to float
data2 = data2.astype(float)

# concatenate data1 (gene list) and data2 (columns of SP and IP) side by side into one dataframe
data = pd.concat([data1, data2], axis=1)
print("data (input data to gsea function): \n", data)
# gene_sets: already in gsea module
gene_sets='GO_Biological_Process_2023'

# run gsea
gs_res = gp.gsea(data=data, # gene expression data table, Pandas Dataframe, gct file
                 gene_sets=gene_sets, # Enrichr library name or .gmt gene sets file or dict of gene sets
                 cls=class_vector, # a list of .cls file
                 # set permutation_type to phenotype if samples >=15
                 permutation_type='phenotype', # choose from "phenotype"(permutate sample labels) or"gene.set" (permutate gene labels)
                 permutation_num=1000, # reduce number to speed up test  # default is 1000
                 outdir=None,  # do not write output to disk  # results output directory
                 method='signal_to_noise',
                 threads=4, seed= 7) # number of threads to use, random seed

"""
:param method:  The method used to calculate a correlation or ranking. Default: 'log2_ratio_of_classes'.
                       Others methods are:
        
                       1. 'signal_to_noise'
        
                          You must have at least three samples for each phenotype to use this metric.
                          The larger the signal-to-noise ratio, the larger the differences of the means
                          (scaled by the standard deviations); that is, the more distinct
                          the gene expression is in each phenotype and the more the gene acts as a “class marker.”
        
                       2. 't_test'
        
                          Uses the difference of means scaled by the standard deviation and number of samples.
                          Note: You must have at least three samples for each phenotype to use this metric.
                          The larger the tTest ratio, the more distinct the gene expression is in each phenotype
                          and the more the gene acts as a “class marker.”
        
                       3. 'ratio_of_classes' (also referred to as fold change).
        
                          Uses the ratio of class means to calculate fold change for natural scale data.
        
                       4. 'diff_of_classes'
        
        
                          Uses the difference of class means to calculate fold change for nature scale data
        
        
                       5. 'log2_ratio_of_classes'
        
                          Uses the log2 ratio of class means to calculate fold change for natural scale data.
                          This is the recommended statistic for calculating fold change for log scale data.
"""


terms = gs_res.res2d.Term  # terms, i.e. enriched pathways
axs = gs_res.plot(terms[:10], show_ranking=False, legend_kws={'loc': (1.05, 0)}, )  # plot top 10 enriched pathways (terms)


"""
GSEA return value:
return: Return a GSEA obj. All results store to a dictionary, obj.results,
             where contains::
    
                 | {
                 |  term: gene set name,
                 |  es: enrichment score,
                 |  nes: normalized enrichment score,
                 |  pval:  Nominal p-value (from the null distribution of the gene set,
                 |  fdr: FDR qvalue (adjusted False Discory Rate),
                 |  fwerp: Family wise error rate p-values,
                 |  tag %: Percent of gene set before running enrichment peak (ES),
                 |  gene %: Percent of gene list before running enrichment peak (ES),
                 |  lead_genes: leading edge genes (gene hits before running enrichment peak),
                 |  matched genes: genes matched to the data,
                 | }
"""

gs_res.res2d  # extract info


from gseapy import heatmap
# plotting 10 heatmaps for top 10 Terms in gsea results
# title of plot is Term, vertical axis is Lead_genes, horizontal axis is samples
# top 10 differentially expressed pathways, can see marked difference between SP and IP
# z-score: normalized data with mean of 0 and variance of 1 across rows (0) or columns (1), in this case rows
for i in range(10):
    genes = gs_res.res2d.Lead_genes[i].split(";")
    # Make sure that ``ofname`` is not None, if you want to save your figure to disk
    ax = heatmap(df = gs_res.heatmat.loc[genes], z_score=0, title=terms[i], figsize=(8,3))


gs_res.heatmat.loc[genes]
# values that plot to heatmap


from gseapy import dotplot
# to save your figure, make sure that ``ofname`` is not None
# we make 3 dotplots, with column values of "NOM p-val", "FDR q-val", "FWER p-val"
# column values produce respective color in dotplot
# row values are same, NES normalized enrichment score
ax1 = dotplot(gs_res.res2d[0:20],  # we use only top 20 rows of gsea results
              x="NES",
             column="NOM p-val",
             title='GO_Biological_Process_2023',
             cmap=plt.cm.viridis,
             size=5,
             figsize=(4,10), cutoff=1, top_term = 20)   # focus on 20 top terms (pathways)

ax2 = dotplot(gs_res.res2d[0:20],
              x="NES",
             column="FDR q-val",
             title='GO_Biological_Process_2023',
             cmap=plt.cm.viridis,
             size=5,
             figsize=(4,10), cutoff=1, top_term = 20)

ax3 = dotplot(gs_res.res2d[0:20],
              x="NES",
             column="FWER p-val",
             title='GO_Biological_Process_2023',
             cmap=plt.cm.viridis,
             size=5,
             figsize=(4,10), cutoff=1, top_term = 20)

# tried to draw dotplots for "Tag %" and "Gene %" too, but there was ValueError
# they are not float values and therefore cannot be plotted into dotplots

#ax4 = dotplot(gs_res.res2d,
#             column="Tag %",
#             title='GO_Biological_Process_2023',
#             cmap=plt.cm.viridis,
#             size=5,
#             figsize=(4,10), cutoff=1, top_term = 20)
# ValueError: some value in Tag % could not be typecast to `float`

#ax5 = dotplot(gs_res.res2d,
#             column="Gene %",
#             title='GO_Biological_Process_2023',
#             cmap=plt.cm.viridis,
#             size=5,
#             figsize=(4,10), cutoff=1, top_term = 20)
# ValueError: some value in Gene % could not be typecast to `float`


print(gs_res)
# gs_res is a gseapy.gsea.GSEA object


enr = gp.enrichr(gene_list=data1, gene_sets='GO_Biological_Process_2023')
enr.results
# enrichment analysis results
# inputted all genes in dataset, shows which pathways are highly related to them


# draw barplot based on enrichment results, top 30 terms, ranked based on adjusted p-value
barplot(enr.results, title='GO_Biological_Process_2023', top_term = 30, figsize = (4,10), color="blue")
plt.show()





enr_map = gp.plot.enrichment_map(gs_res.res2d, column= 'Adjusted P-value', cutoff = 0.2, top_term = 25)
# by changing parameters cutoff and top_term we can change the type and number of elements in the resulting network figure
# cutoff: nodes with 'column' value < cut-off are shown
# -> works only for {'Adjusted P-value', 'P-value', 'NOM p-val', 'FDR q-val'}
# top_term: number of top enriched terms are selected as nodes
print(enr_map)
# return tuple of dataframe (nodes, edges)

# first part is nodes, second part is edges
# nodes is very similar to the gseapy results
# edges contains edge list
# edges show mapping between src_idx to targ_idx, with src_name and targ_name
# and jaccard_coef, overlap_coef and overlap_genes

# Jaccard_coef = intersection /(divided by) union - how strongly related the two pathways are


# networkx docs - https://networkx.org/

import networkx as nx
# import networkx libary
# use it to create network graphs


nodes, edges = enr_map
# retrieve nodes and edges from enrichment map and print them
print("nodes: \n", nodes)
print("\n\n\n\n\nedges: \n", edges)


# build graph
G = nx.from_pandas_edgelist(edges, source='src_idx', target='targ_idx', edge_attr=['jaccard_coef', 'overlap_coef', 'overlap_genes'])
# Returns a graph from 'edges' (edge list)
print(G)


print("G.nodes: \n", G.nodes)  # nodes in graph - nodes in edges
print("G.edges: \n", G.edges)  # list of tuples - two nodes in one edge
nodes_crop = nodes.iloc[list(G.nodes)]  # nodes in graph
print("nodes_crop: \n", nodes_crop)


# We will draw lots of networkx graphs, varying the layout ('pos'). All other code except 'pos' will remain the same

# matplotlib.pyplot.subplots: Create a figure and a set of subplots.
# returns figure and axis (fig, ax)
fig, ax = plt.subplots(figsize=(20, 20))  # figure size
# init node cooridnates
pos=nx.layout.spiral_layout(G)   # lay out nodes in a spiral fashion

# draw nodes
# cmap: Colormap for mapping intensities of nodes
# node size corresponds to the percentage of gene overlap in a certain term of interest
# color of node corresponds to the significance of enriched terms (NES - normalized enrichment score)
# hits: overlapping genes, hits_ratio: ratio of how many genes overlap
nx.draw_networkx_nodes(G, pos=pos, cmap=plt.cm.RdYlBu, node_color=list(nodes_crop.NES), node_size=list(nodes_crop.Hits_ratio *1000))
# draw node label - Terms (pathways) are labels of nodes
nx.draw_networkx_labels(G, pos=pos, labels=nodes_crop.Term.to_dict())

# draw edge
# edge weight is jaccard coefficients of edges
edge_weight = nx.get_edge_attributes(G, 'jaccard_coef').values()

# edge size coresponds to the number of genes that overlap between two connected nodes
# line width of edge signifies jaccard coefficient, degree of overlap between two nodes
nx.draw_networkx_edges(G, pos=pos, width=list(map(lambda x: x*10, edge_weight)), edge_color='#CDDBD4')
plt.show()


# matplotlib.pyplot.subplots: Create a figure and a set of subplots.
# returns figure and axis (fig, ax)
fig, ax = plt.subplots(figsize=(20, 20))  # figure size
# init node cooridnates
pos=nx.layout.circular_layout(G)   # lay out nodes in a circular fashion
# draw node
nx.draw_networkx_nodes(G, pos=pos, cmap=plt.cm.RdYlBu, node_color=list(nodes_crop.NES), node_size=list(nodes_crop.Hits_ratio *1000))
# draw node label
nx.draw_networkx_labels(G, pos=pos, labels=nodes_crop.Term.to_dict())
# draw edge
edge_weight = nx.get_edge_attributes(G, 'jaccard_coef').values()
nx.draw_networkx_edges(G, pos=pos, width=list(map(lambda x: x*10, edge_weight)), edge_color='#CDDBD4')
plt.show()


# Kamada kawai layout reference:
# Tomihisa Kamada and Satoru Kawai. An Algorithm for Drawing General Undirected Graphs. Information Processing Letters 31:7-15, 1989.

# matplotlib.pyplot.subplots: Create a figure and a set of subplots.
# returns figure and axis (fig, ax)
fig, ax = plt.subplots(figsize=(20, 20))  # figure size
# init node cooridnates
pos=nx.layout.kamada_kawai_layout(G)   # Position nodes using Kamada-Kawai path-length cost-function
# draw node
nx.draw_networkx_nodes(G, pos=pos, cmap=plt.cm.RdYlBu, node_color=list(nodes_crop.NES), node_size=list(nodes_crop.Hits_ratio *1000))
# draw node label
nx.draw_networkx_labels(G, pos=pos, labels=nodes_crop.Term.to_dict())
# draw edge
edge_weight = nx.get_edge_attributes(G, 'jaccard_coef').values()
nx.draw_networkx_edges(G, pos=pos, width=list(map(lambda x: x*10, edge_weight)), edge_color='#CDDBD4')
plt.show()


# matplotlib.pyplot.subplots: Create a figure and a set of subplots.
# returns figure and axis (fig, ax)
fig, ax = plt.subplots(figsize=(20, 20))  # figure size
# init node cooridnates
pos=nx.layout.random_layout(G)   # lay out nodes in a random fashion
# draw node
nx.draw_networkx_nodes(G, pos=pos, cmap=plt.cm.RdYlBu, node_color=list(nodes_crop.NES), node_size=list(nodes_crop.Hits_ratio *1000))
# draw node label
nx.draw_networkx_labels(G, pos=pos, labels=nodes_crop.Term.to_dict())
# draw edge
edge_weight = nx.get_edge_attributes(G, 'jaccard_coef').values()
nx.draw_networkx_edges(G, pos=pos, width=list(map(lambda x: x*10, edge_weight)), edge_color='#CDDBD4')
plt.show()


# matplotlib.pyplot.subplots: Create a figure and a set of subplots.
# returns figure and axis (fig, ax)
fig, ax = plt.subplots(figsize=(20, 20))  # figure size
# init node cooridnates
pos=nx.layout.spring_layout(G)
# Position nodes using Fruchterman-Reingold force-directed algorithm.
# The algorithm simulates a force-directed representation of the network treating edges as springs holding nodes close, 
# while treating nodes as repelling objects, sometimes called an anti-gravity force. 
# Simulation continues until the positions are close to an equilibrium.
    
# draw node
nx.draw_networkx_nodes(G, pos=pos, cmap=plt.cm.RdYlBu, node_color=list(nodes_crop.NES), node_size=list(nodes_crop.Hits_ratio *1000))
# draw node label
nx.draw_networkx_labels(G, pos=pos, labels=nodes_crop.Term.to_dict())
# draw edge
edge_weight = nx.get_edge_attributes(G, 'jaccard_coef').values()
nx.draw_networkx_edges(G, pos=pos, width=list(map(lambda x: x*10, edge_weight)), edge_color='#CDDBD4')
plt.show()


# matplotlib.pyplot.subplots: Create a figure and a set of subplots.
# returns figure and axis (fig, ax)
fig, ax = plt.subplots(figsize=(20, 20))  # figure size
# init node cooridnates
pos=nx.layout.spectral_layout(G)
# Position nodes using the eigenvectors of the graph Laplacian.
# Using the unnormalized Laplacian, the layout shows possible clusters of nodes which are an approximation of the ratio cut.
# If dim is the number of dimensions then the positions are the entries of the dim eigenvectors corresponding to the ascending eigenvalues starting from the second one.

# draw node
nx.draw_networkx_nodes(G, pos=pos, cmap=plt.cm.RdYlBu, node_color=list(nodes_crop.NES), node_size=list(nodes_crop.Hits_ratio *1000))
# draw node label
nx.draw_networkx_labels(G, pos=pos, labels=nodes_crop.Term.to_dict())
# draw edge
edge_weight = nx.get_edge_attributes(G, 'jaccard_coef').values()
nx.draw_networkx_edges(G, pos=pos, width=list(map(lambda x: x*10, edge_weight)), edge_color='#CDDBD4')
plt.show()





# gseaplot and gseaplot2 seems to be for Prerank module, not GSEA, therefore of no use
# ringplot is deprecated, use dotplot instead





# What can we do with the gseapy results? Below shows the results
# https://matplotlib.org/stable/plot_types/index.html#plot-types


gs_res.res2d


enr.results


gs_res.res2d[["Term", "NOM p-val"]].transpose()


# Make vertical boxplot
fig1, ax1 = plt.subplots()  # make subplot
ax1.set_title('Vertical Boxplot')  # set title
ax1.boxplot(gs_res.res2d[["ES", "NES", "NOM p-val", "FDR q-val", "FWER p-val"]])  # input data, make boxplot
ax1.set_xticklabels(["ES", "NES", "NOM p-val", "FDR q-val", "FWER p-val"])  # x axis labels


# make horizontal boxplot
fig2, ax2 = plt.subplots()
ax2.set_title('Horizontal Boxplot')
ax2.boxplot(gs_res.res2d[["ES", "NES", "NOM p-val", "FDR q-val", "FWER p-val"]], vert=False)
ax2.set_yticklabels(["ES", "NES", "NOM p-val", "FDR q-val", "FWER p-val"])  # y axis labels


print(enr_map)
# we want to make barplot of hits ratio


# we will make another enrichment map, consisting of top 50 terms at cutoff 0.2, so bigger than the previous one we made
# we will make barplots based on this map

big_enr_map = gp.plot.enrichment_map(gs_res.res2d, column= 'Adjusted P-value', cutoff = 0.2, top_term = 50)
# by changing parameters cutoff and top_term we can change the type and number of elements in the resulting network figure
# cutoff: nodes with 'column' value < cut-off are shown
# -> works only for {'Adjusted P-value', 'P-value', 'NOM p-val', 'FDR q-val'}
# top_term: number of top enriched terms are selected as nodes
print(big_enr_map)
# return tuple of dataframe (nodes, edges)

# first part is nodes, second part is edges
# nodes is very similar to the gseapy results
# edges contains edge list
# edges show mapping between src_idx to targ_idx, with src_name and targ_name
# and jaccard_coef, overlap_coef and overlap_genes

# Jaccard_coef = intersection /(divided by) union - how strongly related the two pathways are


big_nodes, big_edges = big_enr_map
# retrieve nodes and edges from enrichment map and print them


# make another dataframe consisting of columns 'term' and 'hits ratio'
print(big_nodes[["Term", "Hits_ratio"]])
sort_big_nodes = big_nodes[["Term", "Hits_ratio"]].sort_values(by=["Hits_ratio"])
print(sort_big_nodes)


# make boxplot of hits ratio of pathways (terms)

height = sort_big_nodes["Hits_ratio"]
bars = sort_big_nodes["Term"]
y_pos = np.arange(len(bars))

plt.figure(figsize=(20,20))
# Create horizontal bars
plt.barh(y_pos, height)
 
# Create names on the x-axis
plt.yticks(y_pos, bars)

plt.title("Hits ratio of pathways (terms)")

# Show graphic
plt.show()

# Hits: overlapped gene names
# which pathways are most related to our input gene list


# same thing for making another boxplot for normalized enrichment scores (NES)

print(big_nodes[["Term", "NES"]])
nes_big_nodes = big_nodes[["Term", "NES"]].sort_values(by=["NES"])
print(nes_big_nodes)


# create dataset
height = nes_big_nodes["NES"]
bars = nes_big_nodes["Term"]
y_pos = np.arange(len(bars))

plt.figure(figsize=(20,20))
# Create horizontal bars
plt.barh(y_pos, height)
 
# Create names on the x-axis
plt.yticks(y_pos, bars)

plt.title("Normalized enrichment score of pathways (terms)")

# Show graphic
plt.show()

# NES: how enriched the pathway is
# which pathways are most enriched in our dataset










