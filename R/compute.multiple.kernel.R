# compute and returns the multiple kernel
"multiple.kernel" = function( x ) {
    
    # set the parameters
    kernel.type = list()
    kernel.type[1] = list("poly")
    kernel.params = list()
    kernel.params[1] = list(0)
    
    # compute the kernel
    kernels = compute.multiple.kernel(kernel.type,x,x,kernel.params)
    
    # compute some parameters from the kernels
    N = dim(kernels)[1]
    KK = dim(kernels)[3]
    sigma = seq(2,1,-0.25)
    
    # compute and sort Diff
    Diff = dist2(x)^2
    Diff_sort = t(apply(Diff,MARGIN=2,FUN=sort))
    
    # compute the combined kernels
    m = dim(Diff)[1]
    n = dim(Diff)[2]
    allk = seq(10,30,2)
    t = 1
    new_kernels = array(0,c(dim(kernels)[1],dim(kernels)[2],(dim(kernels)[3]+(length(allk)*length(sigma)))))
    new_kernels[,,1:dim(kernels)[3]] = kernels
    kernels = new_kernels
    for (l in 1:length(allk)) {
        if(allk[l]<nrow(x)) {
            TT = t(apply(Diff_sort[,2:(allk[l]+1)],MARGIN=2,FUN=mean))
            Sig = apply(array(0,c(nrow(TT),ncol(TT))),MARGIN=2,FUN=function(x) {x=n})
            Sig = Sig + apply(array(0,c(nrow(TT),ncol(TT))),MARGIN=1,FUN=function(x) {x=n})
            Sig = Sig / 2
            Sig = Sig * which(Sig > .Machine$double.eps,arr.ind=TRUE) + .Machine$double.eps
            for (j in 1:length(sigma)) {
                W = dnorm(Diff,0,sigma[j]*Sig)
                kernels[,,KK+t] = (W + t(W)) / 2
                t = t + 1
            }
        }
    }
    
    # compute D_Kernels
    D_Kernels = kernels
    for (i in 1:dim(kernels)[3]) {
        K = kernels[,,i]
        k = 1/sqrt(diag(K)+1)
        G = K * (k %*% t(k))
        G1 = apply(array(0,c(length(diag(G)),length(diag(G)))),MARGIN=2,FUN=function(x) {x=diag(G)})
        G2 = apply(array(0,c(length(diag(G)),length(diag(G)))),MARGIN=1,FUN=function(x) {x=diag(G)})
        D_Kernels[,,i] = (G1 + G2 - 2*G)/2
        D_Kernels[,,i] = D_Kernels[,,i] - diag(diag(D_Kernels[,,i]))
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
    kernels = array(0,c(n1,n2,k.count))
    
    # compute the kernels
    iteration.kernels = 1
    for (i in 1:count.kernel.type) {
        for (j in 1:length(kernel.params[[i]])) {
            single.kernel.type = kernel.type[[i]]
            single.kernel.parameters = kernel.params[[i]][j]
            if(flag==3) {
                kernels[,,iteration.kernels] = compute.kernel(single.kernel.type,x1,NA,single.kernel.parameters)
            }
            else {
                kernels[,,iteration.kernels] = compute.kernel(single.kernel.type,x1,x2,single.kernel.parameters)
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
    dist = (array(1,c(n2,1)) %*% t(apply(t(x^2),MARGIN=2,FUN=sum))) + 
           (array(1,c(n1,1)) %*% apply(t(c^2),MARGIN=2,FUN=sum)) - 
           - 2 * (x%*%t(c))
    
    return(dist)
    
}
