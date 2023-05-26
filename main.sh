#!/usr/sh
#######################################
#main.sh
#usage:
#sh main.sh \
#	-s human / mouse \
#	-m matrix_file \
#	-l list_file \
#	-c list_column \
#	-a annotation_table
#	-k n_cluster for k-means mathod
#ex.): \
#sh main.sh \
#-s mouse \
#-m Figure.5E-F_Figure.5-1_testicular_gene_expr_table.csv \
#-l Figure.5C_summary_table.csv \
#-c mSymbol \
#-a annotation.csv \
#-k 5
#########################################
#Defining option
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

wget -P input http://purl.obolibrary.org/obo/go/go-basic.obo
wget -P input http://release.geneontology.org/2013-11-01/annotations/goa_${s}.gaf.gz

python Hub_Explorer.py ${m} ${l} ${c} ${a} ${k}

exit=0
