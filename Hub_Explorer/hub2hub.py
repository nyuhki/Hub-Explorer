"""Module Description

    Copyright (c) 2023, Noguchi Yuki

    This code is free software; you can redistribute it and/or modify it
    with following citation:
    Noguchi, Y., Maruoka, M., Onodera, Y., Kosako, H., Suzuki., J. <Title>. <Journal>. <Year>; <Issue>: <Pages>. <DOI>
    @author:  Noguchi Yuki
    @contact: noguchi.yuuki.47u@hotmail.co.jp
"""

import numpy as np
import pandas as pd
from multiprocessing import Pool, get_context, cpu_count

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
