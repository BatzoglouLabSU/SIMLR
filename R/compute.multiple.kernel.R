# compute and returns the multiple kernel
"multiple.kernel" = function( x, cores.ratio = 1 ) {
    
    # set the parameters
    kernel.type = list()
    kernel.type[1] = list("poly")
    kernel.params = list()
    kernel.params[1] = list(0)
    
    # compute some parameters from the kernels
    N = dim(x)[1]
    KK = 0
    sigma = seq(2,1,-0.25)
    
    # compute and sort Diff
    Diff = dist2(x)^2
    Diff_sort = t(apply(Diff,MARGIN=2,FUN=sort))
    
    # compute the combined kernels
    m = dim(Diff)[1]
    n = dim(Diff)[2]
    allk = seq(10,30,2)
    
    # setup a parallelized estimation of the kernels
    cores = as.integer(cores.ratio * (detectCores() - 1))
    if (cores < 1 || is.na(cores) || is.null(cores)) {
        cores = 1
    }

    cl = makeCluster(cores)
    
    clusterEvalQ(cl, {library(Matrix)})
    
    D_Kernels = list()
    D_Kernels = unlist(parLapply(cl,1:length(allk),fun=function(l,x_fun=x,Diff_sort_fun=Diff_sort,allk_fun=allk,
                                                                Diff_fun=Diff,sigma_fun=sigma,KK_fun=KK) {
        if(allk_fun[l]<(nrow(x_fun)-1)) {
            TT = apply(Diff_sort_fun[,2:(allk_fun[l]+1)],MARGIN=1,FUN=mean) + .Machine$double.eps
            TT = matrix(data = TT, nrow = length(TT), ncol = 1)
            Sig = apply(array(0,c(nrow(TT),ncol(TT))),MARGIN=1,FUN=function(x) {x=TT[,1]})
            Sig = Sig + t(Sig)
            Sig = Sig / 2
            Sig_valid = array(0,c(nrow(Sig),ncol(Sig)))
            Sig_valid[which(Sig > .Machine$double.eps,arr.ind=TRUE)] = 1
            Sig = Sig * Sig_valid + .Machine$double.eps
            for (j in 1:length(sigma_fun)) {
                W = dnorm(Diff_fun,0,sigma_fun[j]*Sig)
                D_Kernels[[KK_fun+l+j]] = Matrix((W + t(W)) / 2, sparse=TRUE, doDiag=FALSE)
            }
            return(D_Kernels)
        }
    }))
    
    stopCluster(cl)
    
    # compute D_Kernels
    for (i in 1:length(D_Kernels)) {
        K = D_Kernels[[i]]
        k = 1/sqrt(diag(K)+1)
        G = K * (k %*% t(k))
        G1 = apply(array(0,c(length(diag(G)),length(diag(G)))),MARGIN=2,FUN=function(x) {x=diag(G)})
        G2 = t(G1)
        D_Kernels_tmp = (G1 + G2 - 2*G)/2
        D_Kernels_tmp = D_Kernels_tmp - diag(diag(D_Kernels_tmp))
        D_Kernels[[i]] = Matrix(D_Kernels_tmp, sparse=TRUE, doDiag=FALSE)
    }
    
    return(D_Kernels)
    
}

# compute the multiple kernel
"compute.multiple.kernel" = function( kernel.type, x1, x2 = NA, kernel.params = NA ) {
    
    # set the parameters for x1
    n1 = nrow(x1)
    d1 = ncol(x1)
    
    # set the parameters for x2
    if(any(is.na(x2))) {
        n2 = n1
        d2 = d1
        flag = 0
    }
    else {
        n2 = nrow(x2)
        d2 = ncol(x2)
        if(d1!=d2) {
            stop("Error in the data.")
        }
        flag = 1
    }
    
    # set the count.kernel.type
    count.kernel.type = length(kernel.params)
    if(length(kernel.type)<1 || length(kernel.type)!=count.kernel.type) {
        stop("Error in the parameters.")
    }
    
    # count the parameters
    k.count = 0
    for (i in 1:count.kernel.type) {
        k.count = k.count + length(kernel.params[[i]])
    }
    
    # set the structure to save the kernels
    kernels = list()
    
    # compute the kernels
    iteration.kernels = 1
    for (i in 1:count.kernel.type) {
        for (j in 1:length(kernel.params[[i]])) {
            single.kernel.type = kernel.type[[i]]
            single.kernel.parameters = kernel.params[[i]][j]
            if(flag==3) {
                kernels[[iteration.kernels]] = compute.kernel(single.kernel.type,x1,NA,single.kernel.parameters)
            }
            else {
                kernels[[iteration.kernels]] = compute.kernel(single.kernel.type,x1,x2,single.kernel.parameters)
            }
        }
        iteration.kernels = iteration.kernels + 1
    }
    
    return(kernels)
    
}

# compute the single kernel
"compute.kernel" = function( kernel.type, x1, x2 = NA, kernel.params = NA ) {
    
    # set the parameters for x1
    n1 = nrow(x1)
    d1 = ncol(x1)
    
    # set the parameters for x2
    if(any(is.na(x2))) {
        n2 = n1
        d2 = d1
        flag = 0
    }
    else {
        n2 = nrow(x2)
        d2 = ncol(x2)
        if(d1!=d2) {
            stop("Error in the data.")
        }
        flag = 1
    }
    
    # consider any kernel type
    if(kernel.type=="linear") {
        K = x1 %*% t(x2)
    }
    else if(kernel.type=="poly") {
        K = (x1 %*% t(x2))^kernel.params
    }
    else if(kernel.type=="rbf") {
        if(flag==0) {
            P = apply((x1*x1),MARGIN=1,FUN=sum)
            P1 = t(P)
            P1 = apply(array(0,c(nrow(P1),ncol(P1))),MARGIN=2,FUN=function(x) {x=P1})
            P2 = P
            P2 = apply(array(0,c(nrow(P2),ncol(P2))),MARGIN=1,FUN=function(x) {x=P2})
        }
        else {
            P1 = apply((x1*x1),MARGIN=1,FUN=sum)
            P1 = apply(array(0,c(nrow(P1),ncol(P1))),MARGIN=2,FUN=function(x) {x=P1})
            P2 = t(apply((x2*x2),MARGIN=1,FUN=sum))
            P2 = apply(array(0,c(nrow(P2),ncol(P2))),MARGIN=1,FUN=function(x) {x=P2})
        }
        K = exp(-(P1 + P2 - 2*x1%*%t(x2))%/%(2*kernel.params^2))
    }
    
    return(K)
    
}

# compute the single kernel
"dist2" = function( x, c = NA ) {
    
    # set the parameters for x
    if(is.na(c)) {
        c = x
    }
    
    # compute the dimension
    n1 = nrow(x)
    d1 = ncol(x)
    n2 = nrow(c)
    d2 = ncol(c)
    if(d1!=d2) {
        stop("Data dimension does not match dimension of centres.")
    }
    
    # compute the distance
    dist = t(rep(1,n2) %*% t(apply(t(x^2),MARGIN=2,FUN=sum))) + 
           (rep(1,n1) %*% t(apply(t(c^2),MARGIN=2,FUN=sum))) - 
           2 * (x%*%t(c))
    
    return(dist)
    
}
