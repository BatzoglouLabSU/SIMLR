# compute and returns the multiple kernel
"multiple.kernel.cimlr" = function( x, cores.ratio = 1 ) {
    
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
    Diff = dist2.cimlr(x)
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

# compute the single kernel
"dist2.cimlr" = function( x, c = NA ) {
    
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
    
    dist[which(dist<0,arr.ind=TRUE)] = 0
    
    return(dist)
    
}

# normalizes a symmetric kernel
"dn.cimlr" = function( w, type ) {
    
    # compute the sum of any column
    w = w * dim(w)[1]
    D = apply(abs(w),MARGIN=1,FUN=sum)
    
    # type "ave" returns D^-1*W
    if(type=="ave") {
        D = 1 / D
        D_temp = matrix(0,nrow=length(D),ncol=length(D))
        D_temp[cbind(1:length(D),1:length(D))] = D
        D = D_temp
        wn = D %*% w
    }
    else {
        stop("Invalid type!")
    }
    
    return(wn)
    
}

# umkl function
"umkl.cimlr" = function( D, beta = NA ) {
    
    # set some parameters
    if(is.na(beta)) {
        beta = 1 / length(D)
    }
    tol = 1e-4
    u = 150
    logU = log(u)
    
    # compute Hbeta
    res_hbeta = Hbeta(D, beta)
    H = res_hbeta$H
    thisP = res_hbeta$P
    
    betamin = -Inf
    betamax = Inf
    # evaluate whether the perplexity is within tolerance
    Hdiff = H - logU
    tries = 0
    while (abs(Hdiff) > tol && tries < 30) {
        #if not, increase or decrease precision
        if (Hdiff > 0) {
            betamin = beta
            if(abs(betamax)==Inf) {
                beta = beta * 2
            }
            else {
                beta = (beta + betamax) / 2
            }
        }
        else {
            betamax = beta
            if(abs(betamin)==Inf) {
                beta = beta / 2
            }
            else {
                beta = (beta + betamin) / 2
            }
        }
        # compute the new values
        res_hbeta = Hbeta(D, beta)
        H = res_hbeta$H
        thisP = res_hbeta$P
        Hdiff = H - logU
        tries = tries + 1
    }
    
    return(thisP)
    
}
