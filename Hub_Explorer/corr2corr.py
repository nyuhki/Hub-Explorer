"""Module Description

    Copyright (c) 2023, Noguchi Yuki

    This code is free software; you can redistribute it and/or modify it
    with following citation:
    Noguchi, Y., Maruoka, M., Onodera, Y., Kosako, H., Suzuki., J. <Title>. <Journal>. <Year>; <Issue>: <Pages>. <DOI>
    @author:  Noguchi Yuki
    @contact: noguchi.yuuki.47u@hotmail.co.jp
"""

import pandas as pd
from scipy.stats import spearmanr
from multiprocessing import Pool, get_context, cpu_count

def corr2corr(expr_matrix, corr_list):
	for trial in range(len(corr_list)):
		gene_1 = corr_list[trial]
		opposite_list = expr_matrix.columns.to_list()
		opposite_list.remove(gene_1)
		for gene_2 in opposite_list:
			correlation, pvalue = spearmanr(expr_matrix[gene_1], expr_matrix[gene_2])
			f=open("corr2corr.txt", "a", encoding='utf-8', newline='\n')
			f.write(f'{gene_1},{gene_2},{correlation} \n')
			f.close()
		print(f'{1+trial}/{len(corr_list)}: {gene_1} has been done')

def multi_corr2corr(expr_matrix_dir, corr_list_dir, columns_name):
	#for intact output from Scanpy
	try:
		expr_matrix = pd.read_csv(expr_matrix_dir, index_col=0).drop(["groups"],axis=1)
	except KeyError:
		expr_matrix = pd.read_csv(expr_matrix_dir, index_col=0)
	corr_list = pd.read_csv(corr_list_dir)[columns_name].tolist()
	with get_context("fork").Pool(processes=cpu_count()-1) as pl:
		pl.map(corr2corr(expr_matrix, corr_list), corr_list, 16)
		pl.close()
		pl.join()
        
