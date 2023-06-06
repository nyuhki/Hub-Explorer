# Hub-Explorer
Hub-Explorer aims to unveil the regulatory hub of candidate genes or proteins in the input gene expression matrix likely from development- and temporal biochemical reaction- derived (single cell)RNA-seq data by following three steps:
1. Generation of a hidden gene co-expression (GCN) network:
According to the input expression matrix, the Spearman correlation coefficient between targets and every gene was calculated by the Python pandas.DataFrame.corr function81 in a multiprocessing manner to extract highly-correlated pairs (threshold: spearman ρ≧0.8) as a hidden gene co-expression network (GCN) behind the targets.
2. Hub components identification and back-annotation to target genes:
Each GCN was applied to gene ontology (GO) analysis iteratively using goatools74 to extract significant GO terms (bh_fdr<0.05) as regulatory hub components, and these significant components were back-annotated to the targets.
3. Similarity-based K-means clustering for core regulatory hub identification:
The Jaccard Similarity Index was calculated as described below for evaluating similarity among targets according to the significant hub components:
Jaccard Similarity Index (target A,target B)=  (|target A ∩ target B|)/(|target A ∪ target B|)
Finally, the Python sklearn.cluster.KMeans function86 was utilized to cluster the targets according to the Jaccard Similarity Index, and the overlapped hub components were extracted as the core regulatory hub. As a benchmark for the K-means method, a dendrogram-based clustering was conducted and labeled on the Jaccard Similarity Index-based clustermap by the Python seaborn.clustermap function.

Dependency): \
-Python 3.8.12 \
-goatools v1.2.4 \
-Matplotlib v 3.7.0 \
-Numpy v1.21.2 \
-Pandas v1.5.3 \
-Seaborn v0.12.1 \
-Scikit-learn v0.0.post1 \
-Scipy v1.8.0 

Referece): \
1.Klopfenstein, D.V., Zhang, L., Pedersen, B.S., Ramírez, F., Warwick Vesztrocy, A., Naldi, A., Mungall, C.J., Yunes, J.M., Botvinnik, O.B., Weigel, M., et al. GOATOOLS: A Python library for Gene Ontology analyses. Scientific Reports. 2018; 8: 10872. 10.1038/s41598-018-28948-z \
2.Hunter, J.D. Matplotlib: A 2D Graphics Environment. Computing in Science & Engineering. 2007; 9: 90-95. 10.1109/MCSE.2007.55. \
3.Walt, S.V., Colbert, S.C., and Varoquaux, G. The NumPy Array: A Structure for Efficient Numerical Computation. Computing in Science & Engineering. 2011; 13: 22-30. 10.1109/MCSE.2011.37. \
4.McKinney, W. pandas: a Foundational Python Library for Data Analysis and Statistics. 2011. https://dlr.de/sc/portaldata/15/resources/dokumente/pyhpc2011/submissions/pyhpc2011_submission_9.pdf \
5.Waskom, M.L. Seaborn: Statistical Data Visualization. J. Open Source Softw. 2020; 6: 3021. 10.21105/joss.03021 \
6.Pedregosa, F., Varoquaux, G., Gramfort, A., Michel, V., Thirion, B., Grisel, O., Blondel, M., Louppe, G., Prettenhofer, P., Weiss, R., et al. Scikit-learn: Machine Learning in Python. Journal of Machine Learning Research. 2011; 12: 2825−2830. http://jmlr.org/papers/v12/pedregosa11a.html \
7.Virtanen, P., Gommers, R., Oliphant, T.E., Haberland, M., Reddy, T., Cournapeau, D., Burovski, E., Peterson, P., Weckesser, W., Bright, J., et al. SciPy 1.0: fundamental algorithms for scientific computing in Python. Nature Methods. 2019; 17: 261 - 272. 10.1038/s41592-019-0686-2




