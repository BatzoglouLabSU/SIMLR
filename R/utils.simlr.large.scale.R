# function to compute fast PCA
"fast.pca" = function( x, k ) {
    
    # performs svds
    res_svds = svds(x,k)
    U = res_svds$u
    S = res_svds$d
    
    # perform fast PCA
    tmp = array(0,c(length(S),length(S)))
    diag(tmp) = S
    res = U %*% tmp
    
    return(res)
    
}

# function to perform k-nearest neighbour search and distances
"knn.search" = function( xx, n ) {
    
    # performs k-nearest neighbour search
    ind = randomProjectionTreeSearch(xx,n,distance_method="Euclidean")
    
    # compute the distances from the results
    indices = neighborsToVectors(ind)
    val = distance(i=indices$i,j=indices$j,x=xx,distance_method="Euclidean")
    
    res = list()
    res$val = t(matrix(abs(val),nrow=n,ncol=dim(xx)[2]))
    res$ind = ind + 1
    
    return(res)
    
}

# compute and returns the multiple kernel for large scale data
"multiple.kernel_large_scale" = function( val, ind, kk = 20 ) {
    
    # compute some parameters from the kernels
    sigma = seq(2,1,-0.25)
    
    # compute the combined kernels
    allk = seq(from=ceiling(kk/2),to=ceiling(kk*1.5),by=ceiling(kk/10))
    
    D_Kernels = list()
    val_fun=val
    ind_fun=ind
    allk_fun=allk
    sigma_fun=sigma
    KK_fun=0
    for(l in 1:length(allk)) {
        if(allk_fun[l]<(ncol(val_fun))) {
            TT = apply(val_fun[,1:allk_fun[l]],MARGIN=1,FUN=mean) + .Machine$double.eps
            TT0 = apply(array(0,c(ncol(val_fun))),MARGIN=1,FUN=function(x) {x=TT})
            TT0 = (TT0 + TT[ind_fun])*0.5
            for (j in 1:length(sigma_fun)) {
                temp = dnorm(val_fun,0,sigma_fun[j]*TT0)
                temptemp = temp[,1]
                temp = (apply(array(0,c(ncol(val_fun))),MARGIN=1,FUN=function(x) {x=temptemp}) + temptemp[ind_fun]) * 0.5 - temp
                KK_fun = KK_fun + 1
                D_Kernels[[KK_fun]] = temp + .Machine$double.eps
            }
        }
    }
    
    return(D_Kernels)
    
}

# compute the L2 distance for large scale datasets
"L2_distance_large_scale" <- function( a, b ) {
    
    I = matrix(rep(1:dim(b)[1],dim(b)[2]),nrow=dim(b)[1],ncol=dim(b)[2])
    temp = apply((F[as.vector(I),]-F[b,])^2,MARGIN=2,FUN=sum)
    
    d = matrix(temp,nrow=dim(b)[1],ncol=dim(b)[1])
    
    return(d)
    
}
