"""Module Description

    Copyright (c) 2023, Noguchi Yuki

    This code is free software; you can redistribute it and/or modify it
    with following citation:
    Noguchi, Y., Maruoka, M., Onodera, Y., Kosako, H., Suzuki., J. <Title>. <Journal>. <Year>; <Issue>: <Pages>. <DOI>
    @author:  Noguchi Yuki
    @contact: noguchi.yuuki.47u@hotmail.co.jp

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

"""

import sys
import glob
from Hub_Explorer.corr2corr import *
from Hub_Explorer.corr2hub_components import *
from Hub_Explorer.hub2hub import *
from Hub_Explorer.hub_explorer import *

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
	hub_explorer(data_dir = glob.glob("./Module_Summary/Summary_module*.csv"))
