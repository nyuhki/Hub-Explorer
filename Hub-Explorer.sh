#!/usr/sh
###########################################################################
##Hub-Explorer.sh
#
#    Copyright (c) 2023, Noguchi Yuki
#
#    This code is free software; you can redistribute it and/or modify it
#    with following citation:
#    Noguchi, Y., Maruoka, M., Onodera, Y., Kosako, H., Suzuki., J. <Title>. <Journal>. <Year>; <Issue>: <Pages>. <DOI>
#    @author:  Noguchi Yuki
#    @contact: noguchi.yuuki.47u@hotmail.co.jp
#
####################
##   DESCRIPTION  ##
####################
##USAGE (INPUT):
#sh Hub-Explorer.sh \
#	-s human / mouse \
#	-m matrix_file \
#	-l list_file \
#	-c list_column \
#	-a annotation_table
#	-k n_cluster for k-means mathod
#
##OUTPUT:
#Hidden regulatory hub behinds of interest protein / gene list
#
##ADVANTAGE:
#1.Infer the function of unannotated genes and novel isoforms such as TE-isoform.
#2.Dissect the temporal gene dynamics with a development-based dataset such as a spermatogenesis scRNA-seq.
#
##PROCESS:
#	1. Generation of a hidden gene co-expression (GCN) network
#		According to the input expression matrix, the Spearman correlation coefficient was calculated
#		by the "pandas.DataFrame.corr" function between each protein and every gene.
#		The highly-correlated pair (threshold: spearman ρ = 0.8) was extracted as the hidden gene
#		co-expression network (GCN) behind each input protein.
#	2. Back-annotation of the regulatory hub components
#		Each GCN was applied to gene ontology (GO) analysis iteratively using goatools,
#		and significant GO terms (bh_fdr < 0.05) were extracted as the regulatory hub components
#		and back-annotated to the input protein.
#	3. Extraction of the regulatory hub based on K-mean clustering
#		According to the hub components, the Jaccard Similarity Index was calculated among input proteins
#		following the formula:
#				Jaccard Similarity Index (Protein A,Protein B)=  (|Protein A ∩ Protein B|)/(|Protein A ∪ Protein B|)
#	Finally, the "sklearn.cluster.KMeans" function was applied to the similarity matrix to group input proteins
#	and extract the overlapped hub components as the regulatory hub. As a benchmark for the k-mean method,
#	the dendrogram-based clustering was performed and overlapped on the heatmap.
#
##REFERENCES
#1.	Klopfenstein, D.V., Zhang, L., Pedersen, B.S., Ramírez, F., Warwick Vesztrocy, A., Naldi, A., Mungall, C.J., Yunes, J.M., Botvinnik, O.B., Weigel, M., et al. GOATOOLS: A Python library for Gene Ontology analyses. Scientific Reports. 2018; 8: 10872. 10.1038/s41598-018-28948-z
###########################################################################

function usage {
    cat <<EOM
Usage: $(basename "$0") [OPTION]...
    -h          Display help
    -s VALUE    Species (Homo Sapience:human, Mus Musculus:mouse)
    -m VALUE    matrix_file.csv
    -l VALUE    list_file.csv
    -c VALUE    column_name
    -a VALUE    annotation_table.csv
    -k VALUE    n_cluster
EOM
    exit 2
}

#Definition
while getopts ":s:m:l:c:a:k:h" optKey; do
    case "$optKey" in
        s)
          s=$(echo "${OPTARG}")
          ;;
		m)
          m=$(echo "${OPTARG}")
          ;;
		l)
          l=$(echo "${OPTARG}")
          ;;
		c)
          c=$(echo "${OPTARG}")
          ;;
		a)
          a=$(echo "${OPTARG}")
          ;;
		k)
          k=$(echo "${OPTARG}")
          ;;
        '-h'|'--help'|* )
          usage
          ;;
    esac
done

rm ./input/*.obo ./input/*.gaf.gz
wget -P ./input http://purl.obolibrary.org/obo/go/go-basic.obo
if [ "${s}" = "human" ]; then
	wget -P ./input http://http.ebi.ac.uk/pub/databases/GO/goa/HUMAN/goa_human.gaf.gz
else
	wget -P ./input http://http.ebi.ac.uk/pub/databases/GO/goa/MOUSE/goa_mouse.gaf.gz
fi

python __main__.py ${s} ${m} ${l} ${c} ${a} ${k}

mkdir result
mv *.txt ./result/
mv *.csv ./result/
mv *.eps ./result/
mv ./GO_result/ ./result/
mv ./Module_GO/ ./result/
mv ./Module_Summary ./result/

exit=0

