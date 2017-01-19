# perform the SIMLR clustering algorithm for large scale datasets
"SIMLR_Large_Scale" <- function( X, c, k = 10, kk = 100, if.impute = FALSE, normalize = FALSE, cores.ratio = 1 ) {
    
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
    
    # start the clock to measure the execution time
    ptm = proc.time()
    
    # set some parameters
    NITER = 5
    num = ncol(X)
    r = -1
    beta = 0.8
    X = t(X)
    
    cat("Performing fast PCA.\n")
    fast.pca_res = fast.pca(X,k=kk)
    fast.pca_res = t(fast.pca_res)
    
    cat("Performing k-nearest neighbour search.\n")
    nearest_neighbour_res = knn.search(fast.pca_res,k*3)
    val = nearest_neighbour_res$val
    ind = nearest_neighbour_res$ind
    
    cat("Computing the multiple Kernels.\n")
    
    # compute the kernels
    D_Kernels = multiple.kernel_large_scale(val,ind,k)
    
    # set up some parameters
    alphaK = 1 / rep(length(D_Kernels),length(D_Kernels))
    distX = array(0,c(dim(D_Kernels[[1]])[1],dim(D_Kernels[[1]])[2]))
    for (i in 1:length(D_Kernels)) {
        distX = distX + D_Kernels[[i]]
    }
    distX = distX / length(D_Kernels)
    
    # sort distX for rows
    res = apply(distX,MARGIN=1,FUN=function(x) return(sort(x,index.return = TRUE)))
    distX1 = array(0,c(nrow(distX),ncol(distX)))
    idx = array(0,c(nrow(distX),ncol(distX)))
    for(i in 1:nrow(distX)) {
        distX1[i,] = res[[i]]$x
        idx[i,] = res[[i]]$ix
    }
    
    di = distX1[,2:(k+2)]
    rr = 0.5 * (k * di[,k+1] - apply(di[,1:k],MARGIN=1,FUN=sum))
    
    if(r<=0) {
        r = mean(rr)
    }
    lambda = max(mean(rr),0)
    S0 = max(max(distX)) - distX
    
    cat("Performing the iterative procedure ",NITER," times.\n")
    
    # compute dn
    S0 = dn_large_scale(S0,'ave')
    
    S0_sparse = sparseMatrix(i=as.vector(matrix(rep(1:nrow(ind),ncol(ind)))),j=as.vector(ind),x=as.vector(S0))
    eig_res = eigs_sym(S0_sparse,c)
    F_eig = eig_res$vectors
    eig_res = eig_res$values
    
    # perform the iterative procedure NITER times
    for(iter in 1:NITER) {
        
        cat("Iteration: ",iter,"\n")
        
        distf = L2_distance_large_scale(F_eig,ind)
        b = idx[,2:dim(idx)[2]]
        a = apply(array(0,c(num,ncol(b))),MARGIN=2,FUN=function(x){ x = 1:num })
        inda = cbind(as.vector(a),as.vector(b))
        ad = (distX[inda]+lambda*distf[inda])/2/r
        dim(ad) = c(num,ncol(b))
        
        # call the c function for the optimization
        c_input = -t(ad)
        c_output = t(ad)
        ad = t(.Call("projsplx_R",c_input,c_output))
        
        ad = dn(S0,'ave')
        S0 = beta * S0 + (1 - beta) * ad
        
        S0_sparse = sparseMatrix(i=as.vector(matrix(rep(1:nrow(ind),ncol(ind)))),j=as.vector(ind),x=as.vector(S0))
        eig_res = eigs_sym(S0_sparse,c)
        F_eig = eig_res$vectors
        eig_res = eig_res$values
        
        DD = vector()
        for (i in 1:length(D_Kernels)) {
            temp = (.Machine$double.eps+D_Kernels[[i]]) * (S0+.Machine$double.eps)
            DD[i] = mean(apply(temp,MARGIN=2,FUN=sum))
        }
        alphaK0 = umkl(DD)
        alphaK0 = alphaK0 / sum(alphaK0)
        alphaK = (1-beta) * alphaK + beta * alphaK0
        alphaK = alphaK / sum(alphaK)
        lambda = 1.5 * lambda
        r = r / 1.01
        
        # compute Kbeta
        distX = D_Kernels[[1]] * alphaK[1]
        for (i in 2:length(D_Kernels)) {
            distX = distX + as.matrix(D_Kernels[[i]]) * alphaK[i]
        }
        
    }
    
    # compute the execution time
    execution.time = proc.time() - ptm
    
    cat("Performing Kmeans.\n")
    y = kmeans(F_last,c,nstart=200)
    
    ydata = NULL
    
    # create the structure with the results
    results = list()
    results[["y"]] = Y
    results[["S0"]] = S0
    results[["F"]] = F_eig
    results[["ydata"]] = ydata
    results[["alphaK"]] = alphaK
    results[["val"]] = val
    results[["ind"]] = ind
    results[["execution.time"]] = execution.time
    
    return(results)
    
}
