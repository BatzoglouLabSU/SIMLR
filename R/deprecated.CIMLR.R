#' CIMLR_Estimate_Number_of_Clusters
#' 
#' @title This function has been removed and it is now provided
#' by the CIMLR package, please refer to https://github.com/danro9685/CIMLR
#'
#' @examples    
#' \dontrun{    
#' CIMLR_Estimate_Number_of_Clusters(GliomasReduced$in_X,   
#'    NUMC = 2:5,   
#'    cores.ratio = 0)  
#' }
#' @param all_data is a list of multi-omic data each of which is an (m x n) data matrix of measurements of cancer patients  
#' @param NUMC vector of number of clusters to be considered    
#' @param cores.ratio ratio of the number of cores to be used when computing the multi-kernel   
#'  
#' @export CIMLR_Estimate_Number_of_Clusters
#' @return a list of 2 elements: K1 and K2 with an estimation of the best clusters (the lower   
#' values the better) as discussed in the original paper of SIMLR
#'
"CIMLR_Estimate_Number_of_Clusters" = function(
    all_data,
    NUMC = 2:5,
    cores.ratio = 1) {
    if (!requireNamespace("CIMLR", quietly = TRUE)) {
        stop("Package \"CIMLR\" needed for this function to work. Please install it from https://github.com/danro9685/CIMLR",
        call. = FALSE)
    }
    .Defunct("CIMLR_Estimate_Number_of_Clusters", package = "CIMLR", 
        msg = "This function has been removed and it is now provided by the CIMLR package, please refer to https://github.com/danro9685/CIMLR")
}

#' perform the CIMLR clustering algorithm   
#'  
#' @title This function has been removed and it is now provided
#' by the CIMLR package, please refer to https://github.com/danro9685/CIMLR 
#'  
#' @examples    
#' \dontrun{
#' CIMLR(X = GliomasReduced$in_X, c = 3, cores.ratio = 0)   
#' }    
#' @param X a list of multi-omic data each of which is an (m x n) data matrix of measurements of cancer patients    
#' @param c number of clusters to be estimated over X   
#' @param no.dim number of dimensions   
#' @param k tuning parameter    
#' @param cores.ratio ratio of the number of cores to be used when computing the multi-kernel   
#'  
#' @return clusters the patients based on CIMLR and their similarities  
#' @export CIMLR        
#'  
"CIMLR" <- function( X, c, no.dim = NA, k = 10, cores.ratio = 1 ) {
    if (!requireNamespace("CIMLR", quietly = TRUE)) {
        stop("Package \"CIMLR\" needed for this function to work. Please install it from https://github.com/danro9685/CIMLR",
        call. = FALSE)
    }
    .Defunct("CIMLR", package = "CIMLR", 
        msg = "This function has been removed and it is now provided by the CIMLR package, please refer to https://github.com/danro9685/CIMLR")
}
