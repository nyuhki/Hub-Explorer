"""Module Description

    Copyright (c) 2023, Noguchi Yuki

    This code is free software; you can redistribute it and/or modify it
    with following citation:
    Noguchi, Y., Maruoka, M., Onodera, Y., Kosako, H., Suzuki., J. <Title>. <Journal>. <Year>; <Issue>: <Pages>. <DOI>
    @author:  Noguchi Yuki
    @contact: noguchi.yuuki.47u@hotmail.co.jp
"""

import os
import shutil
import pandas as pd
import numpy as np
from goatools import obo_parser
import Bio.UniProt.GOA as GOA
import gzip
from goatools.go_enrichment import GOEnrichmentStudy

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
  with gzip.open(f'./input/goa_{species}.gaf.gz', 'rt') as gaf_fp:
    funcs = {}  # Initialise the dictionary of functions
    # Iterate on each function using Bio.UniProt.GOA library.
    for entry in GOA.gafiterator(gaf_fp):
        uniprot_id = entry.pop('DB_Object_ID')
        funcs[uniprot_id] = entry
  pop = funcs.keys()
  assoc = {}
  for x in funcs:
    if x not in assoc:
      assoc[x] = set()
    assoc[x].add(str(funcs[x]['GO_ID']))
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
