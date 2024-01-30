#!/usr/sh
###########################################################################
#
#Current Version is Deposited to `https://github.com/SuzukiLab-icems/Hub-Explorer`
#Future update will be conducted in this Repository.
#
#Hub-Explorer_v1.0.0/Hub-Explorer.sh
#
#	 Copyright (c) 2024, Noguchi Yuki (Jun Suzuki lab)
#	 This software is released under the MIT License, see LICENSE (https://opensource.org/license/mit/).
#    @citation: Noguchi, Y., Onodera, Y., Maruoka, M., Miyamoto, T., Kosako, H., Suzuki., J. 2024. In vivo CRISPR screening directly targeting testicular cells. Cell Genomics.
#    @author:  Noguchi Yuki
#    @contact: nyuhki21@gmail.com,jsuzuki@icems.kyoto-u.ac.jp
#
##REFERENCE
#1.	Klopfenstein, D.V., Zhang, L., Pedersen, B.S., Ramírez, F., Warwick Vesztrocy, A., Naldi, A., Mungall, C.J., Yunes, J.M., Botvinnik, O.B., Weigel, M., et al. GOATOOLS: A Python library for Gene Ontology analyses. Scientific Reports. 2018; 8: 10872. 10.1038/s41598-018-28948-z
#
##CMD for generating figures
#Figure.5E-5G: sh ./Hub-Explorer_v1.0.0/Hub-Explorer.sh -i IO -m mtx_from_Figure.4B-C_testicular_gene_expr_table_for_Figure.5.csv -l list_modified_from_Figure.5C_summary_table_of_spermatid_specific_RD3_interactors.csv -a annotation_file_from_scanpy_var.csv -o manual -k 4 -t Rd3
#Figure.S5D-S5E: sh ./Hub-Explorer_v1.0.0/Hub-Explorer.sh -i IO -m mtx_from_Figure.4B-C_testicular_gene_expr_table_for_Figure.5.csv -l list_from_Figure.S5C.csv -a annotation_file_from_scanpy_var.csv -o manual -k 7 -t None
#
#
##For test run
#You can run Hub-Explorer by typing 'sh Hub-Explorer_v1.0.0/Hub-Explorer.sh -i demo_data -m demo_matrix.csv -l demo_list.csv -a annotation_file.csv -o manual -k 2 -t None' in 'Hub-Explorer' directory. Please type pwd and confirm your current directory is 'Hub-Explorer.'
#
#[Directory Architecture]
#Hub-Explorer
#	|
#	|-Hub-Explorer_v1.0.0
#	|	|-source codes..
#	|
#	|-demo_data
#		|-annotation_file.csv
#		|-demo_list.csv
#		|-demo_matrix.csv
#		|-goa_mouse.gaf.gz
#		|-go-basic.obo
#
#If you want to change the directory, pls fix line.107 according to your preference.
###########################################################################
