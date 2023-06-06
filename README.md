# Hub-Explorer
Hub-Explorer aims to unveil the regulatory hub of candidate genes / proteins in the specific manner like a spermatogenesis, retinogenesis and cancer clonal evolution 
1. Generation of a hidden gene co-expression (GCN) network:
According to the input expression matrix, the Spearman correlation coefficient between targets and every gene was calculated by the Python pandas.DataFrame.corr function81 in a multiprocessing manner to extract highly-correlated pairs (threshold: spearman ρ≧0.8) as a hidden gene co-expression network (GCN) behind the targets.
2. Hub components identification and back-annotation to target genes:
Each GCN was applied to gene ontology (GO) analysis iteratively using goatools74 to extract significant GO terms (bh_fdr<0.05) as regulatory hub components, and these significant components were back-annotated to the targets.
3. Similarity-based K-means clustering for core regulatory hub identification:
The Jaccard Similarity Index was calculated as described below for evaluating similarity among targets according to the significant hub components:
Jaccard Similarity Index (target A,target B)=  (|target A ∩ target B|)/(|target A ∪ target B|)
Finally, the Python sklearn.cluster.KMeans function86 was utilized to cluster the targets according to the Jaccard Similarity Index, and the overlapped hub components were extracted as the core regulatory hub. As a benchmark for the K-means method, a dendrogram-based clustering was conducted and labeled on the Jaccard Similarity Index-based clustermap by the Python seaborn.clustermap function.
