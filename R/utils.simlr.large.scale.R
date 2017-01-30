# perform fast pca
"fast.pca" = function( X, K ) {

    X = t(X)
    tmp_val = as.vector(colSums(X)/nrow(X))
    X = X - t(apply(array(0,dim(X)),MARGIN=1,FUN=function(x) {x=tmp_val}))
    res = fast.rsvd(X,K)
    U = res$U
    S = res$S
    K = min(dim(S)[2],K)
    diag_val = sqrt(diag(S[1:K,1:K]))
    diag_mat = array(0,c(length(diag_val),length(diag_val)))
    diag(diag_mat) = diag_val
    X = U[,1:K]%*%diag_mat
    normalization_val = sqrt(rowSums(X*X))
    X = X / apply(array(0,c(length(normalization_val),K)),MARGIN=2,FUN=function(x) {x=normalization_val})

    return(X)

}

# perform fast rsvd
"fast.rsvd" = function( A, K ) {
    
    M = dim(A)[1]
    N = dim(A)[2]
    P = min(2*K,N)
    X = matrix(rnorm(N*P),nrow=N,ncol=P)
    Y = A%*%X
    W1 = orth(Y)
    B = t(W1)%*%A
    res = svd(B,nu=min(dim(B)),nv=min(dim(B)))
    W2 = res$u
    tmp_S = res$d
    S = array(0,c(length(tmp_S),length(tmp_S)))
    diag(S) = tmp_S
    V = res$v
    U = W1%*%W2
    K = min(K,dim(U)[2])
    U = U[,1:K]
    S = S[1:K,1:K]
    V = V[,1:K]

    return(list(U=U,S=S,V=V))

}

# compute and returns the multiple kernel for large scale data
"multiple.kernel_large_scale" = function( val, ind, kk = 20 ) {
    
    # compute some parameters from the kernels
    sigma = seq(2,1,-0.25)
    
    # compute the combined kernels
    allk = seq(from=ceiling(kk/2),to=ceiling(kk*1.5),by=ceiling(kk/10))
    
    D_Kernels = list()
    val_fun=val*val
    ind_fun=ind
    allk_fun=allk
    sigma_fun=sigma
    KK_fun=0
    for(l in 1:length(allk)) {
        if(allk_fun[l]<(ncol(val_fun))) {
            TT = apply(val_fun[,1:allk_fun[l]],MARGIN=1,FUN=mean) + .Machine$double.eps
            TT0 = apply(array(0,c(ncol(val_fun))),MARGIN=1,FUN=function(x) {x=TT})
            TT0 = (TT0 + matrix(TT[ind_fun],nrow=dim(TT0)[1],ncol=dim(TT0)[2]))*0.5
            for (j in 1:length(sigma_fun)) {
                temp = dnorm(val_fun,0,sigma_fun[j]*TT0)
                temptemp = temp[,1]
                temp = (apply(array(0,c(ncol(val_fun))),MARGIN=1,FUN=function(x) {x=temptemp}) + matrix(temptemp[ind_fun],nrow=dim(TT0)[1],ncol=dim(TT0)[2])) * 0.5 - temp
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
    temp = rowSums((a[as.vector(I),]-a[as.vector(b),])^2)
    
    d = matrix(temp,nrow=dim(b)[1],ncol=dim(b)[2])
    
    return(d)
    
}

# normalizes a symmetric kernel for large scale
"dn_large_scale" = function( w, type ) {
    
    # compute the sum of any column
    D = apply(abs(w),MARGIN=1,FUN=sum)
    
    # type "ave" returns D^-1*W
    if(type=="ave") {
        D = 1 / D
        D_temp = matrix(0,nrow=length(D),ncol=length(D))
        D_temp[cbind(1:length(D),1:length(D))] = D
        D = D_temp
        wn = D %*% w
    }
    # type "gph" returns D^-1/2*W*D^-1/2
    else if(type=="gph") {
        D = 1 / sqrt(D)
        D_temp = matrix(0,nrow=length(D),ncol=length(D))
        D_temp[cbind(1:length(D),1:length(D))] = D
        D = D_temp
        wn = D %*% (w %*% D)
    }
    else {
        stop("Invalid type!")
    }
    
    return(wn)
    
}
