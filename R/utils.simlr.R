# compute the eigenvalues and eigenvectors
"eig1" <- function( A, c = NA, isMax = NA, isSym = NA ) {
    
    # set the needed parameters
    if(is.na(c)) {
        c = dim(A)[1]
    }
    if(c>dim(A)[1]) {
        c = dim(A)[1]
    }
    if(is.na(isMax)) {
        isMax = 1
    }
    if(is.na(isSym)) {
        isSym = 1
    }
    
    # compute the eigenvalues and eigenvectors of A
    if(isSym==1) {
        eigen_A = eigen(A,symmetric=TRUE)
    }
    else {
        eigen_A = eigen(A)
    }
    v = eigen_A$vectors
    d = eigen_A$values
    
    # sort the eigenvectors
    if(isMax == 0) {
        eigen_A_sorted = sort(d,index.return=TRUE)
    }
    else {
        eigen_A_sorted = sort(d,decreasing=TRUE,index.return=TRUE)
    }
    d1 = eigen_A_sorted$x
    idx = eigen_A_sorted$ix
    idx1 = idx[1:c]
    
    # compute the results
    eigval = d[idx1]
    eigvec = Re(v[,idx1])
    eigval_full = d[idx]
    
    return(list(eigval=eigval,eigvec=eigvec,eigval_full=eigval_full))
    
}

# compute the L2 distance
"L2_distance_1" <- function( a, b ) {
    
    if(dim(a)[1] == 1) {
        a = rbind(a,rep(0,dim(a)[2]))
        b = rbind(b,rep(0,dim(b)[2]))
    }
    
    aa = apply(a*a,MARGIN=2,FUN=sum)
    bb = apply(b*b,MARGIN=2,FUN=sum)
    ab = t(a) %*% b
    d1 = apply(array(0,c(length(t(aa)),length(bb))),MARGIN=2,FUN=function(x){ x = t(aa) })
    d2 = t(apply(array(0,c(length(t(bb)),length(aa))),MARGIN=2,FUN=function(x){ x = t(bb) }))
    d = d1 + d2 - 2 * ab
    d = Re(d)
    d = matrix(mapply(d,FUN=function(x) { return(max(max(x),0)) }),nrow=nrow(d),ncol=ncol(d))
    d_eye = array(1,dim(d))
    diag(d_eye) = 0
    d = d * d_eye
    
    return(d)
    
}

# umkl function
"umkl" = function( D, beta = NA ) {
    
    # set some parameters
    if(is.na(beta)) {
        beta = 1 / length(D)
    }
    tol = 1e-4
    u = 20
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
                beta = beta * 2
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

"Hbeta" = function( D, beta ) {
    
    D = (D - min(D)) / (max(D) - min(D) + .Machine$double.eps)
    P = exp(-D * beta)
    sumP = sum(P)
    H = log(sumP) + beta * sum(D * P) / sumP
    P = P / sumP
    
    return(list(H=H,P=P))
    
}
