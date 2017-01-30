# perform the SIMLR clustering algorithm for large scale datasets
"SIMLR_Large_Scale" <- function( X, c, k = 10, kk = 100, if.impute = FALSE, normalize = FALSE ) {
    
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
    
    cat("Performing fast PCA.\n")
    
    fast.pca_res = fast.pca(X,kk)
    
    cat("Performing k-nearest neighbour search.\n")
    
    a_annoy = new(AnnoyEuclidean,dim(fast.pca_res)[2])
    n_annoy = dim(fast.pca_res)[1]
    for (i_annoy in seq(n_annoy)) {
        v_annoy = as.vector(fast.pca_res[i_annoy,])
        a_annoy$addItem(i_annoy-1,v_annoy)
    }
    a_annoy$build(200)
    val = array(0,c(dim(fast.pca_res)[1],k*2))
    ind = array(0,c(dim(fast.pca_res)[1],k*2))
    for(j_annoy in 1:dim(val)[1]) {
        ind[j_annoy,] = a_annoy$getNNsByItem(j_annoy-1,k*2)+1
        for(val_annoy in 1:dim(val)[2]) {
            val[j_annoy,val_annoy] = a_annoy$getDistance(j_annoy-1,ind[j_annoy,val_annoy]-1)
        }
    }
    
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
    
    di = distX[,2:(k+2)]
    rr = 0.5 * (k * di[,k+1] - apply(di[,1:k],MARGIN=1,FUN=sum))
    
    if(r<=0) {
        r = mean(rr)
    }
    lambda = max(mean(rr),0)
    S0 = max(max(distX)) - distX
    
    # compute dn
    S0 = dn_large_scale(S0,'ave')
    
    S0_sparse = sparseMatrix(i=as.vector(matrix(rep(1:nrow(ind),ncol(ind)))),j=as.vector(ind),x=as.vector(S0),dims=c(nrow(ind),nrow(ind)))
    eig_res = eigs_sym(S0_sparse+t(S0_sparse),c,which="LM")
    F_eig = Re(eig_res$vectors)
    eig_res = Re(eig_res$values)/max(Re(eig_res$values))
    eig_res = (1-beta)*eig_res / (1-beta*eig_res^2)
    tmp_eig = array(0,c(length(eig_res),length(eig_res)))
    diag(tmp_eig) = eig_res
    F_eig = F_eig %*% tmp_eig
    F_eig = dn_large_scale(F_eig,'ave')
    F_eig_initial = F_eig
    
    cat("Performing the iterative procedure ",NITER," times.\n")
    
    # perform the iterative procedure NITER times
    for(iter in 1:NITER) {
        
        cat("Iteration: ",iter,"\n")
        
        distf = L2_distance_large_scale(F_eig,ind)
        ad = (distX+lambda*distf)/2/r
        
        # call the c function for the optimization
        c_input = -t(ad)
        c_output = t(ad)
        ad = t(.Call("projsplx_R",c_input,c_output))
        
        S0 = (1 - beta) * S0 + beta * ad
        
        S0_sparse = sparseMatrix(i=as.vector(matrix(rep(1:nrow(ind),ncol(ind)))),j=as.vector(ind),x=as.vector(S0),dims=c(nrow(ind),nrow(ind)))
        eig_res = eigs_sym(S0_sparse+t(S0_sparse),c,which="LM")
        F_eig = Re(eig_res$vectors)
        eig_res = Re(eig_res$values)/max(Re(eig_res$values))
        eig_res = (1-beta)*eig_res / (1-beta*eig_res^2)
        tmp_eig = array(0,c(length(eig_res),length(eig_res)))
        diag(tmp_eig) = eig_res
        F_eig = F_eig %*% tmp_eig
        F_eig = dn_large_scale(F_eig,'ave')
        F_eig = (1 - beta) * F_eig_initial + beta * F_eig
        F_eig_initial = F_eig
        
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
    
    y = kmeans(F_eig,c,nstart=200)
    
    cat("Performing t-SNE.\n")
    
    I = as.vector(seq(1,ncol(ind)*nrow(ind)+1,ncol(ind)))
    J = as.vector(t(ind))
    V = as.vector(t(S0))
    ydata = SIMLR.Rtsne(I,J,V)$Y
    
    # create the structure with the results
    results = list()
    results[["y"]] = y
    results[["S0"]] = S0
    results[["F"]] = F_eig
    results[["ydata"]] = ydata
    results[["alphaK"]] = alphaK
    results[["val"]] = val
    results[["ind"]] = ind
    results[["execution.time"]] = execution.time
    
    return(results)
    
}
