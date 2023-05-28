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
