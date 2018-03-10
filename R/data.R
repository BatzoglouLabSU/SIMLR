#' @name BuettnerFlorian 
#' @title test dataset for SIMLR 
#' @description example dataset to test SIMLR from the work by Buettner, Florian, et al. 
#' @docType data 
#' @usage data(BuettnerFlorian) 
#' @format gene expression measurements of individual cells 
#' @source Buettner, Florian, et al. "Computational analysis of cell-to-cell heterogeneity in single-cell RNA-sequencing data reveals hidden subpopulations of cells." Nature biotechnology 33.2 (2015): 155-160. 
#' @return list of 6: 
#'		in_X = input dataset as an (m x n) gene expression measurements of individual cells, 
#'  	n_clust = number of clusters (number of distinct true labels), 
#'  	true_labs = ground true of cluster assignments for each of the n_clust clusters, 
#'  	seed = seed used to compute the results for the example, 
#'  	results = result by SIMLR for the inputs defined as described, 
#'  	nmi = normalized mutual information as a measure of the inferred clusters compared to the true labels 
NULL

#' @name ZeiselAmit 
#' @title test dataset for SIMLR large scale 
#' @description example dataset to test SIMLR large scale. This is a reduced version of the dataset from the work by Zeisel, Amit, et al. 
#' @docType data 
#' @usage data(ZeiselAmit) 
#' @format gene expression measurements of individual cells 
#' @source Zeisel, Amit, et al. "Cell types in the mouse cortex and hippocampus revealed by single-cell RNA-seq." Science 347.6226 (2015): 1138-1142. 
#' @return list of 6: 
#'		in_X = input dataset as an (m x n) gene expression measurements of individual cells, 
#'  	n_clust = number of clusters (number of distinct true labels), 
#'  	true_labs = ground true of cluster assignments for each of the n_clust clusters, 
#'  	seed = seed used to compute the results for the example, 
#'  	results = result by SIMLR for the inputs defined as described, 
#'  	nmi = normalized mutual information as a measure of the inferred clusters compared to the true labels 
NULL

#' @name GliomasReduced 
#' @title test dataset for CIMLR 
#' @description example dataset to test CIMLR. This is a reduced version of the dataset from the work by The Cancer Genome Atlas Research Network. 
#' @docType data 
#' @usage data(GliomasReduced) 
#' @format multi-omic data of cancer patients 
#' @source Cancer Genome Atlas Research Network. "Comprehensive, integrative genomic analysis of diffuse lower-grade gliomas." New England Journal of Medicine 372.26 (2015): 2481-2498. 
#' @return list of 1 element: 
#'		in_X = input dataset as a list of 4 (reduced) multi-omic data each of which is an (m x n) measurements of cancer patients 
NULL
