"""
Hub_Explorer.py matrix_file(args[1]) list_file(args[2]) list_column(args[3]) annotation_table(args[4])
ex.):
Hub_Explorer.py \
	mouse \
	Figure.5E-F_Figure.5-1_testicular_gene_expr_table.csv \
	Figure.5C_summary_table(spermatid_specific).csv \
	mSymbol \
	annotation_(scanpy_var).csv

#Input
-species(args[1])			:mouse / human
-matrix_file(args[2])		:Expression data set
-list_file(args[3])			:Table containing target gene / protein name
-list_column(args[4])		:Column name containing target gene / protein name
-annotation_table(args[5])  :Corresponding table between the gene / protein names and uniprot accession IDs
-n_cluster(arg[6])		    :No of cluster for k-mean clustering

#Processing
Calculate the spearman-correlation value between each input target vs. every gene following the input data set (expected pattern: number of targets × 20,000 genes)
Generate a correlation gene list corresponding to each input candidate by extracting highly-correlated pairs (spearman correlation > 0.8)
Perform GO analysis toward the correlated gene list. This GO information is a functional feature of each input candidate.
Annotate the GO information for each input recursively.
Calculate Jaccard Similarity Index among input candidates according to annotated GO information semantically.
Perform k-means clustering to cluster input candidates functionally.
Extract overlapped GO terms in each of the clusters.

#Output
-Functional regulation relationships among targeted gene/protein

#Advantage
-Infer the function of unannotated genes and novel isoforms such as TE-isoform.
-Dissect the temporal gene dynamics with a development-based dataset such as a spermatogenesis scRNA-seq.

#Detailed Description
This process inputs an expression data set, such as scRNA-seq data and a list of target genes/proteins likely from proteome analysis. Next, it calculates the Spearman correlation value between each target and every gene in the data set (number of candidates vs. around 20,000 genes). It then extracts the highly-correlated pair (spearman correlation > 0.8). These correlations generate a correlation gene list specific to each input candidate. The gene ontology (GO) analysis is then performed on the correlation gene list, which provides functional information about each candidate. Finally, this GO information is recursively annotated for each input candidate.
Subsequently, the Jaccard Similarity Index is semantically calculated among the input candidates based on the annotated GO information, measuring the functional similarity among candidates. To further analyze the candidates, k-means clustering is applied, grouping them based on functional similarities. Eventually, the overlapped GO terms within each cluster are extracted.
The output of this process is a functional connectome of candidate proteins/genes constructed based on the identified GO terms. One advantage of this method is its usefulness in inferring the functions of unreviewed genes, such as novel isoforms. Additionally, it can be applied to analyze temporal gene regulation dynamics using development-related datasets, such as spermatogenesis single-cell RNA sequencing (scRNA-seq) data.

#Reference
1.	Klopfenstein, D.V., Zhang, L., Pedersen, B.S., Ramírez, F., Warwick Vesztrocy, A., Naldi, A., Mungall, C.J., Yunes, J.M., Botvinnik, O.B., Weigel, M., et al. GOATOOLS: A Python library for Gene Ontology analyses. Scientific Reports. 2018; 8: 10872. 10.1038/s41598-018-28948-z
"""
import sys
import os
import glob
import re
import collections
import argparse
import pandas as pd
import numpy as np
from pandas.api.types import is_numeric_dtype
import matplotlib.pyplot as plt
import scipy.stats as stats
import seaborn as sns
from scipy.stats import spearmanr
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
import scanpy as sc
from sklearn.cluster import AgglomerativeClustering
from sklearn.cluster import KMeans
from sklearn.metrics import jaccard_score
from goatools import obo_parser
import Bio.UniProt.GOA as GOA
from ftplib import FTP
import gzip
from goatools.go_enrichment import GOEnrichmentStudy
from numba import jit
from multiprocessing import Pool, get_context, cpu_count
from concurrent import futures
from concurrent.futures import ThreadPoolExecutor
import shutil

def corr2corr(expr_matrix, corr_list):
    for gene_1 in corr_list:
        opposite_list = expr_matrix.columns.to_list()
        opposite_list.remove(gene_1)
        for gene_2 in opposite_list:
            correlation, pvalue = spearmanr(expr_matrix[gene_1], expr_matrix[gene_2])
            f=open("corr2corr.txt", "a", encoding='utf-8', newline='\n')
            f.write(f'{gene_1},{gene_2},{correlation} \n')
            f.close()
            print(f'{gene_1} vs. {gene_2}')

def multi_corr2corr(expr_matrix_dir, corr_list_dir, columns_name):
    expr_matrix = pd.read_csv(expr_matrix_dir, index_col=0).drop(["groups"],axis=1)
    corr_list = pd.read_csv(corr_list_dir)[columns_name].tolist()
    with get_context("fork").Pool(processes=cpu_count()-1) as pl:
        pl.map(corr2corr(expr_matrix, corr_list), corr_list, 16)
        pl.close()
        pl.join()

def corr2list():
    corr2corr = pd.read_table('corr2corr.txt', sep=',', usecols=[0,1,2], header=None, names=['gene_1','gene_2','rho'], index_col=0)
    corr2corr = corr2corr[corr2corr["rho"] != "nan "]
    corr2corr["rho"] = corr2corr["rho"].astype("float")
    corr2corr_list = corr2corr[corr2corr["rho"] > 0.8]
    return corr2corr_list

def matrix2go(corr2corr_list, annotation_file, obo_file, methods, species):
  if os.path.isdir("./GO_result/"):
    shutil.rmtree("./GO_result/")
    print("GO_result dir exists.")
  os.mkdir("./GO_result/")
  go = obo_parser.GODag(obo_file)
  with gzip.open(f'./input/goa_{species}.gaf.gz', 'rt') as mouse_gaf_fp:
    mouse_funcs = {}  # Initialise the dictionary of functions
    # Iterate on each function using Bio.UniProt.GOA library.
    for entry in GOA.gafiterator(mouse_gaf_fp):
        uniprot_id = entry.pop('DB_Object_ID')
        mouse_funcs[uniprot_id] = entry
  pop = mouse_funcs.keys()
  assoc = {}
  for x in mouse_funcs:
    if x not in assoc:
      assoc[x] = set()
    assoc[x].add(str(mouse_funcs[x]['GO_ID']))
  #Prep for study file
  annotation = \
    pd.read_csv(annotation_file)\
    .rename(columns={"Unnamed: 0":"Gene symbol", "gene_ids":"Ensembl_ID"})\
    .loc[:,["Gene symbol","Uniprot ID"]]
  corr2corr = corr2corr_list.reset_index().rename(columns={"gene_2":"Gene symbol"})
  for gene in np.unique(corr2corr["gene_1"].to_list()):
    gene2go = corr2corr[corr2corr["gene_1"] == gene]
    list = gene2go.merge(annotation, how="left", on="Gene symbol")
    list = list.dropna()["Uniprot ID"].tolist()
    study = {}
    for x in list:
      if x not in study:
        study[x] = set()
    study = study.keys()
    #Calculation
    g_fdr = GOEnrichmentStudy(pop, assoc, go,
                              propagate_counts=True,
                              alpha=0.05,
                              methods=methods)
    g_fdr_res = g_fdr.run_study(study)
    g_fdr.wr_xlsx(f'./GO_result/{gene}.xlsx', g_fdr_res)

def raw2valid(raw_data):
    sig_table = pd.DataFrame(columns=["GO", "name", "p_fdr_bh", "study_items", "gene_info"])
    for list in raw_data:
        list = list.replace('.xlsx', '').replace('./GO_result/','')
        tmp = pd.read_excel(f'./GO_result/{list}.xlsx').loc[:,["GO", "name", "p_fdr_bh", "study_items"]]
        tmp = tmp[tmp["p_fdr_bh"] < 0.05]
        tmp["gene_info"] = list
        sig_table = pd.concat([sig_table, tmp])
        print(f'{list}')
    sig_table.dropna().to_csv("valid_te2go_table.csv")

def go2cluster(gene_1):
    def jaccard_similarity(x, y):
        intersection = len(set.intersection(*[set(x), set(y)]))
        union = len(set.union(*[set(x), set(y)]))
        return intersection / float(union)
    sig_table = pd.read_csv("valid_te2go_table.csv")
    x_GO = sig_table[sig_table["gene_info"] == gene_1]["GO"].tolist()
    for gene_2 in np.unique(sig_table["gene_info"].tolist()):
        y_GO = sig_table[sig_table["gene_info"] == gene_2]["GO"].tolist()
        try:
            index = jaccard_similarity(x_GO, y_GO)
            f=open("./Jaccard_Index.csv", "a", encoding='utf-8', newline='\n')
            f.write(f'{gene_1}_vs_{gene_2}, {index: .4f}\n')
            f.close()
            print(f'{gene_1}_vs_{gene_2}')
        except "ZeroDivisionError":
            try:
                f=open("./Jaccard_Index.csv", "a", encoding='utf-8', newline='\n')
                f.write(f'{gene_1}_vs_{gene_2}, ZeroDivisionError\n')
                f.close()
                print(f'{gene_1}_vs_{gene_2}')
            except "TypeError":
                f=open("./Jaccard_Index.csv", "a", encoding='utf-8', newline='\n')
                f.write(f'{gene_1}_vs_{gene_2}, TypeError\n')
                f.close()
                print(f'{gene_1}_vs_{gene_2}')

def multi_go2cluster():
    sig_table = pd.read_csv("valid_te2go_table.csv")
    list = np.unique(sig_table["gene_info"].tolist())
    with get_context("fork").Pool(processes=cpu_count()-1) as pl:
        pl.map(go2cluster, list, 16)
        pl.close()
        pl.join()

def jc2matrix():
    jc = pd.read_csv("./Jaccard_Index.csv", header=None, names=['combination', 'jaccard_index'])
    jc_tmp = jc["combination"].str.split("_vs_", expand=True).set_axis(['gene_1', 'gene_2'], axis='columns')
    jc = jc.reset_index().merge(jc_tmp.reset_index(), how="left", on="index").set_index("combination").drop("index",axis=1)
    jc = jc.reset_index().loc[:,["combination","gene_1","gene_2","jaccard_index"]]
    jc.to_csv("GO_similarity_list.csv")
    jc_matrix = jc.pivot_table(index="gene_1", columns="gene_2").transpose().reset_index().drop("level_0", axis=1).set_index("gene_2")
    jc_matrix.to_csv("GO_similarity_matrix.csv")

def kmeans(corr_matrix, n_clusters):
    Kmodel = KMeans(n_clusters=n_clusters, random_state=0, init='random')
    Kmodel.fit(corr_matrix)
    #Processing
    cluster = Kmodel.labels_.tolist()
    df_annotation = pd.DataFrame(data={'Cluster': cluster}, index=corr_matrix.index).reset_index()
    df_annotation.to_csv("te_cluster_label.csv")
    clustered = df_annotation.sort_values(by="Cluster", ascending=False).merge(corr_matrix.reset_index(), how="left", on="gene_2").set_index("gene_2")
    clustered = df_annotation.rename(columns = {"gene_2":"index"}) \
                             .merge(clustered.drop("Cluster",axis=1).transpose().reset_index(),how="left",on="index") \
                             .sort_values(by="Cluster", ascending=False).drop("Cluster",axis=1).set_index("index")
    annotated_cluster = df_annotation.rename(columns = {"gene_2":"index"}) \
                                     .merge(clustered.reset_index(), how="left", on="index") \
                                     .sort_values(by="Cluster", ascending=False).set_index("index").rename(columns={"index":"Gene"})
    annotated_cluster.to_csv("gene_classified_matrix.csv")
    return annotated_cluster, clustered, df_annotation

def kmeans2GOannotation(df_annotation, n_clusters):
    if os.path.isdir("./Module_GO/"):
        shutil.rmtree("./Module_GO/")
        print("Module_GO dir exists.")
    os.mkdir("./Module_GO/")
    """annotation"""
    go_list = pd.read_csv("valid_te2go_table.csv",index_col=0) \
                .loc[:,["gene_info","GO","name","p_fdr_bh","study_items"]] \
                .rename(columns={"gene_info":"Gene", "study_items":"Uniprot_ID"})
    go2cluster = go_list.merge(df_annotation.rename(columns={"gene_2":"Gene"}), how="left", on="Gene") \
                        .loc[:,["Gene","Cluster","GO","name","p_fdr_bh","Uniprot_ID"]] \
                        .sort_values(by="Cluster", ascending=True)
    go2cluster["Uniprot_ID"] = go2cluster["Uniprot_ID"] + str(",")
    for module in np.unique(go2cluster["Cluster"].tolist()):
        go2cluster_module = go2cluster[go2cluster["Cluster"] == module]
        gene_1 = go2cluster_module["Gene"].tolist()[0]
        x = go2cluster_module[go2cluster_module["Gene"] == gene_1]["GO"].tolist()
        opposite_list = go2cluster_module["Gene"].tolist()
        opposite_list.remove(gene_1)
        shared = set(x)
        """shared_GO"""
        for gene_2 in opposite_list:
            y = go2cluster_module[go2cluster_module["Gene"] == gene_2]["GO"].tolist()
            shared = shared & set(y)
        if len(shared) != 0:
            df_shared = pd.DataFrame(data=shared, columns=["GO"]) \
                          .merge(go_list.loc[:,["GO"]],how="left",on="GO") \
                          .groupby("GO").sum() \
                          .merge(go_list.loc[:,["GO","name"]].drop_duplicates(),how="left",on="GO").loc[:,["GO","name"]]
            df_shared.to_csv(f'./Module_GO/shared_GO_module{module}.csv')
        else:
            print(f'Module{module}: There is no shared GO')
    go2cluster.to_csv("go2cluster.csv")
    return go2cluster

def module_summary(go2cluster, annotated_cluster):
    if os.path.isdir("./Module_Summary/"):
        shutil.rmtree("./Module_Summary/")
        print("Module_Summary dir exists.")
    os.mkdir("./Module_Summary/")
    go2cluster = go2cluster.loc[:,["Cluster","Gene","GO","name","p_fdr_bh","Uniprot_ID"]]
    for module in np.unique(go2cluster["Cluster"].tolist()):
        go2cluster[go2cluster["Cluster"] == module].to_csv(f'./Module_Summary/Summary_module{module}.csv')
    annotated_cluster.reset_index().loc[:,["index","Cluster"]].groupby("Cluster").count().rename(columns={"index":"count"}).to_csv("./Module_Summary/module_count_summary.csv")

def heatmap(annotated_clustered, n_clusters, dendrogram, min, center, max, format_type, show):
    zeileis_colors = np.array(sc.pl.palettes.godsnot_102)
    col_colors = np.array(annotated_clustered["Cluster"]).astype("<U7")
    for color in np.arange(n_clusters):
      col_colors = np.where(col_colors == str(color), zeileis_colors[[color]], col_colors)
    sns.set(style="ticks",
        font_scale=2.5,
        font='arial',
        rc = {'figure.figsize':(15,8), 'axes.linewidth':1.5})
    g = sns.clustermap(
        annotated_clustered.drop("Cluster",axis=1),
        method="ward",
        cmap="GnBu",
        center=center,
        vmin=min, vmax=max,
        col_cluster=dendrogram,
        row_cluster=dendrogram,
        col_colors=col_colors,
        row_colors=col_colors,
        xticklabels=False,
        yticklabels=False)
    g.fig.subplots_adjust(right=0.7)
    g.ax_cbar.set_position((0.8, .2, .03, .4))
    g.ax_cbar.set_title('J.I.', rotation=0)
    plt.tight_layout()
    plt.savefig(f'cluster.{format_type}',format=f'{format_type}',dpi=150)
    if show:
        plt.show()
    plt.close()

def main_clustering2plot(n_clusters, dendrogram, format_type, show):
    """Clustering"""
    jc_matrix = pd.read_csv("./GO_similarity_matrix.csv", index_col=0)
    jc_matrix = jc_matrix.loc[:,jc_matrix.index]
    annotated_cluster, clustered, df_annotation = kmeans(jc_matrix,n_clusters)
    go2cluster = kmeans2GOannotation(df_annotation, n_clusters)
    module_summary(go2cluster, annotated_cluster)
    """Plot"""
    heatmap(annotated_clustered=annotated_cluster, n_clusters=n_clusters, dendrogram=dendrogram, min=0, center=0.375, max=1.0, format_type=format_type, show=show)

def kmeans2jc(module_1):
    def jaccard_similarity(x, y):
        intersection = len(set.intersection(*[set(x), set(y)]))
        union = len(set.union(*[set(x), set(y)]))
        return intersection / float(union)
    go_list = pd.read_csv("kmeans2jc_table.csv",index_col=0)
    x_GO = go_list[go_list["module"] == module_1]["GO"].tolist()
    for module_2 in np.unique(go_list["module"].tolist()):
        y_GO = go_list[go_list["module"] == module_2]["GO"].tolist()
        try:
            index = jaccard_similarity(x_GO, y_GO)
            f=open("./kmeans2similarity.csv", "a", encoding='utf-8', newline='\n')
            f.write(f'{module_1}_vs_{module_2}, {index: .4f}\n')
            f.close()
            print(f'{module_1}_vs_{module_2}')
        except "ZeroDivisionError":
            try:
                f=open("./kmeans2similarity.csv", "a", encoding='utf-8', newline='\n')
                f.write(f'{module_1}_vs_{module_2}, ZeroDivisionError\n')
                f.close()
                print(f'{module_1}_vs_{module_2}')
            except "TypeError":
                f=open("./kmeans2similarity.csv", "a", encoding='utf-8', newline='\n')
                f.write(f'{module_1}_vs_{module_2}, TypeError\n')
                f.close()
                print(f'{module_1}_vs_{module_2}')

def multi_kmeans2jc(n_clusters):
    if os.path.isfile("./kmeans2similarity.csv"):
        os.remove("./kmeans2similarity.csv")
    go_list = pd.DataFrame(columns=["GO","name","module"])
    for module in np.arange(n_clusters):
        try:
            df_module = pd.read_csv(f'./Module_GO/shared_GO_module{module}.csv', index_col=0)
            df_module["module"] = module
            go_list = pd.concat([go_list, df_module])
        except FileNotFoundError:
            print(f'There is no shared module in module{module}.')
    go_list.to_csv("kmeans2jc_table.csv")
    module_1 = np.unique(go_list["module"].tolist())
    with get_context("fork").Pool(processes=cpu_count()-1) as pl:
        pl.map(kmeans2jc, module_1, 16)
        pl.close()
        pl.join()

def heatmap_kmeans2jc(annotated_clustered, n_clusters, dendrogram, min, center, max, format_type, show):
    zeileis_colors = np.array(sc.pl.palettes.godsnot_102)
    col_colors = np.array(annotated_clustered.index).astype("<U7")
    for color in np.array(annotated_clustered.index):
        col_colors = np.where(col_colors == str(color), zeileis_colors[[color]], col_colors)
    sns.set(style="ticks",
        font_scale=2.5,
        font='arial',
        rc = {'figure.figsize':(15,8), 'axes.linewidth':1.5})
    g = sns.clustermap(
        data=annotated_clustered,
        method="ward",
        cmap="mako",
        center=center,
        vmin=min, vmax=max,
        col_cluster=dendrogram,
        row_cluster=dendrogram,
        col_colors=col_colors,
        row_colors=col_colors,
        xticklabels=False,
        yticklabels=False)
    g.fig.subplots_adjust(right=0.7)
    g.ax_cbar.set_position((0.8, .2, .03, .4))
    g.ax_cbar.set_title('J.I.', rotation=0)
    plt.tight_layout()
    plt.savefig(f'GO_cluster.{format_type}',format=f'{format_type}',dpi=150)
    if show:
        plt.show()
    plt.close()

def kmeans2jc_main(n_clusters, format_type, dendrogram, show):
    def kmeans2jc_matrix():
        jc = pd.read_csv("./kmeans2similarity.csv", header=None, names=['combination', 'jaccard_index'])
        jc_tmp = jc["combination"].str.split("_vs_", expand=True).set_axis(['module_A', 'module_B'], axis='columns')
        jc = jc.reset_index().merge(jc_tmp.reset_index(), how="left", on="index").set_index("combination").drop("index",axis=1)
        jc = jc.reset_index().loc[:,["combination","module_A","module_B","jaccard_index"]]
        kmeans2jc_matrix = jc.pivot_table(index="module_A", columns="module_B").transpose().reset_index().drop("level_0", axis=1).set_index("module_B")
        kmeans2jc_matrix.to_csv("kmeans2jc_matrix.csv")
    """Clustering"""
    multi_kmeans2jc(n_clusters)
    kmeans2jc_matrix()
    """Plot"""
    kmeans2jc_matrix = pd.read_csv("kmeans2jc_matrix.csv", index_col=0)
    heatmap_kmeans2jc(annotated_clustered=kmeans2jc_matrix, n_clusters=n_clusters, dendrogram=dendrogram, min=0, center=0.375, max=1.0, format_type=format_type, show=show)

def network_inference(data_dir = glob.glob("./Module_Summary/Summary_module*.csv")):
	def final_heatmap(anno_table, format_type, show):
		sns.set(style="ticks",
			font_scale=2.5,
			font='arial',
			rc = {'figure.figsize':(40,18), 'axes.linewidth':1.5})
		sns.heatmap(
			anno_table.transpose(),
			cmap="Pastel1_r",
			vmin=0, vmax=1,
			cbar=False
			)
		plt.xticks(fontsize=34, rotation=25, ha ='right')
		plt.yticks(fontsize=34)
		plt.tight_layout()
		plt.savefig(f'GO_table.{format_type}',format=f'{format_type}',dpi=150)
		plt.close()
	"""matrix_processing"""
	df = pd.DataFrame(columns=["Cluster","Gene","GO","name"])
	for data in data_dir:
		df_tmp = pd.read_csv(data, usecols=[1,2,3,4])
		df = pd.concat([df,df_tmp])
	df = df.sort_values(by="Cluster", ascending=True).set_index("Cluster")
	df["Count"] = 1
	table = df.pivot_table(index="Gene", columns="name")
	table = table.fillna(0).transpose().reset_index().drop("level_0", axis=1).set_index("name")
	label = df.reset_index().loc[:,["Cluster", "Gene"]]
	label = label.drop_duplicates()
	name_label = df.reset_index().loc[:,["Cluster","name"]]
	name_label = name_label.set_index("Cluster")["name"].drop_duplicates()
	anno_table = label.merge(table.transpose().reset_index(), how="left", on="Gene")
	tmp = pd.DataFrame(columns=anno_table.columns)
	for module in np.unique(anno_table["Cluster"].tolist()):
		df_module = anno_table[anno_table["Cluster"] == module]
		df_module["gene_scale"] = anno_table.sum(axis=1)
		df_module = df_module.sort_values(by="gene_scale", ascending=False).drop("gene_scale",axis=1)
		tmp = pd.concat([tmp, df_module])
	anno_table = tmp
	anno_table = anno_table.drop("Cluster",axis=1).set_index("Gene").transpose()
	anno_table = name_label.reset_index().merge(
				anno_table.reset_index().rename(columns={"index":"name"}),
				how="left",
				on="name").sort_values(by="Cluster", ascending=True)
	frame = pd.DataFrame(columns = anno_table.columns)
	for module in np.unique(anno_table["Cluster"].tolist()):
		df_module = anno_table[anno_table["Cluster"] == module]
		df_module["GO_scale"] = anno_table.sum(axis=1)
		df_module = df_module.sort_values(by="GO_scale", ascending=True).drop("GO_scale",axis=1)
		frame = pd.concat([frame,df_module])
	anno_table = frame.drop("Cluster", axis=1).set_index("name")
	final_heatmap(anno_table, format_type="eps", show=False)

if __name__ == "__main__":
	args = sys.argv
	species = args[1]
	matrix_file = args[2]
	list_file = args[3]
	list_column = args[4]
	annotation_table = args[5]
	n_cluster = int(args[6])
	"""correlation_matrix"""
	try:
		multi_corr2corr(
			expr_matrix_dir = "./input/" + matrix_file,
			corr_list_dir = "./input/" + list_file,
			columns_name = list_column
			)
	except TypeError:
		print(f'corr2corr.txt has been generated.')
	"""matrix2go"""
	corr2corr_list = corr2list()
	try:
		corr2corr_list = corr2corr_list.drop(["Rd3"])
	except:
		corr2corr_list = corr2corr_list
	matrix2go(
        corr2corr_list=corr2corr_list,
        annotation_file = "./input/" + annotation_table,
        obo_file = "./input/go-basic.obo",
        methods = ["fdr_bh"],
        species = species)
	raw2valid(raw_data=glob.glob('./GO_result/*.xlsx'))
	"""go2cluster"""
	multi_go2cluster()
	jc2matrix()
	"""Kmean clustering, JC_heatmap"""
	main_clustering2plot(n_clusters=n_cluster, dendrogram=True, format_type='eps', show=False)
	kmeans2jc_main(n_clusters=n_cluster, format_type='eps', dendrogram=True, show=False)
	"""GO and gene table"""
	network_inference(data_dir = glob.glob("./Module_Summary/Summary_module*.csv"))
