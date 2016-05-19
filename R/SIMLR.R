# perform the SIMLR clustering algorithm
"SIMLR" <- function( X, c , no.dim = NA, k = 10, if.impute = FALSE, normalize = FALSE ) {
	
	# set any required parameter to the defaults
	if(is.na(no.dim)) {
		no.dim = c
	}
	
	# check the if.impute parameter
	if(if.impute == TRUE) {
		X = t(X)
		X_zeros = which(X==0,arr.ind=TRUE)
		if(length(X_zeros)>0) {
			R_zeros = as.vector(X_zeros[,"row"])
			C_zeros = as.vector(X_zeros[,"col"])
			ind = (C_zeros - 1) * nrow(X) + R_zeros
			X[ind] = as.vector(colMeans(X))[C_zeros]
		}
		X = t(X)
	}
	
	# check the normalize parameter
	if(normalize == TRUE) {
		X = t(X)
		X = X - min(as.vector(X))
		X = X / max(as.vector(X))
		C_mean = as.vector(colMeans(X))
		X = apply(X,MARGIN=1,FUN=function(x) return(x-C_mean))
	}
	
}
