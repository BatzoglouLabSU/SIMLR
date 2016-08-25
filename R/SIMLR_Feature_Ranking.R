# A is the similarity matrix by SIMLR
# X is the data of size nxp
"SIMLR_Feature_Ranking" <- function( A, X ) {
    
    yscore = array(NA,c(100,nrow(A)))
    for (i in 1:100) {
        cat(i,"\n")
        index = sample(1:nrow(A))
        index = index[1:round(nrow(A)*0.9)]
        Ai = A[index,index]
        Xi = X[index,]
        res = LaplacianScore(Xi,Ai)
        yscore[i,] = res
    }

    yscore = 1 - yscore
    numerator = (t(yscore) - min(as.vector(yscore)) + .Machine$double.eps)
    denominator = (max(as.vector(yscore)) - min(as.vector(yscore)) + .Machine$double.eps)
    glist = numerator / denominator
    res = aggregateRanks(glist)
    res2 = sort(res$pval,index.return=TRUE)
    
    res = list(pval=res2$x,aggR=res$aggR[res2$ix])
    
    return(res)
    
}

"LaplacianScore" = function( X, W ) {
    
    nSmp = ncol(X)
    nFea = nrow(X)
    
    if(ncol(W) != nFea) {
        stop('W is error')
    }

    D = apply(W,MARGIN=1,FUN=sum)
    L = W

    allone = rep(1,nSmp)

    tmp1 = t(D) %*% X
    
    DPrime = sum(t((t(X)%*%D)))-tmp1*tmp1/sum(diag(D))
    LPrime = sum(t((t(X)%*%L))*X)-tmp1*tmp1/sum(diag(D))
    DPrime[DPrime < 1e-12] = 10000
    Y = LPrime/DPrime
    Y = t(Y)
    
    return(as.vector(Y))
    
}

"aggregateRanks" = function( R, N = NA, method = "RRA", complete = NA, topCutoff = NaN ) {
    
    if(is.null(R)) {
        stop('Input parameter R is missing!')
    }
    
    if(is.na(complete)) {
        complete = 0
    }
    
    if (all(apply(R,MARGIN=2,FUN=max)<=1) && all(apply(R,MARGIN=2,FUN=min)<=1)>0) {
        rmat = R
    }
    else {
        stop('Columns of matrix R can only contain numbers from interval (0,1].')
    }
    
    if(method=="min") {
        aggR = apply(rmat[!is.nan(rmat)],MARGIN=1,FUN=min)
        pval = NA
    }
    else if(method=="median") {
        aggR = apply(rmat[!is.nan(rmat)],MARGIN=1,FUN=median)
        pval = NA
    }
    else if(method=="geom.mean") {
        aggR = apply(rmat[!is.nan(rmat)],MARGIN=1,FUN=function(x) return(exp(mean(log(rmat)))))
        pval = NA
    }
    else if(method=="mean") {
        aggR = apply(rmat[!is.nan(rmat)],MARGIN=1,FUN=mean)
        n = apply(rmat[!is.nan(rmat)],MARGIN=1,FUN=sum)
        pval = dnorm(aggR,mean=0.5,sd=sqrt(1/12/n))
    }
    else if(method=="RRA") {
        aggR = rhoScores(rmat,topCutoff)
        pval = aggR
    }
    else {
        stop('Method should be one of:  "min", "geom.mean", "mean", "median" or "RRA"')
    }
    
    return(list(pval=pval,aggR=aggR))
    
}

"rhoScores" = function( r, topCutoff = NaN ) {
    
    rho = rep(NaN,length(r))
    
    for(rInd in 1:length(r)) {
        r1 = r[rInd]
        if(is.nan(topCutoff)) {
            x = betaScores(r1)
            rho[rInd] = correctBetaPvalues(min(x),sum(x[!is.nan(x)],na.rm=TRUE))
        }
        else{
            r1 = r1[!is.nan(r1)]
            r1[r1==1] = rep(NaN,length(which(r1==1)))
            x = thresholdBetaScore(r=r1,topCutoff=topCutoff)
            rho[rInd] = correctBetaPvalues(min(x),length(r1))
        }
    }
    
    return(rho)
    
}

"betaScores" = function( r ) {
    
    n = sum(r,na.rm=TRUE)
    p = rep(NA,1:length(r))
    r = sort(r)
    p[1:n] = dbeta(r[1:n],1:n,n:1)
    return(p)
    
}

"thresholdBetaScore" = function( r, k = NA, n = NA, sigma = NA ) {
    
    # set the defaults values is needed
    rLen = length(r)
    if(is.na(k)) {
        k = 1:rLen
    }
    if(is.na(n)) {
        n = rLen
    }
    if(is.na(sigma)) {
        sigma = rep(1,rLen)
    }
    
    # check for any error
    if(length(sigma) != n) {
        stop('The length of sigma does not match n!')
    }
    if(length(r) != n) {
        stop('The length of p-values does not match n!')
    }
    if(min(sigma)< 0 || max(sigma) > 1) {
        stop('Elements of sigma are not in the range [0,1]!')
    }
    if(any(!is.nan(r) & r > sigma)) {
        stop('Elements of r must be smaller than elements of sigma!')
    }
    
    x = sort(r[!is.nan(r)])
    sigma = sort(sigma,decreasing=TRUE)
    Beta = rep(NaN,length(k))
    
    for(i in 1:length(k)) {
        
        if(k[i]>n) {
            Beta[i] = 0
        }
        if(k[i]> length(x)) {
            Beta[i] = 1
        }
        if(sigma[n]>= x[k[i]]) {
            Beta[i] = dbeta(x[k[i]],k[i],(n+1-k[i]))
        }
        
        n0 = which(sigma<x[k[i]])
        n0 = n0[1] - 1
        if(n0==0) {
            B = c(1,rep(0,k[i]))
        }
        else if(k[i]>n0) {
            B = c(1,dbeta(x[k[i]],1:n0, n0:1),rep(0,(k[i]-n0)))
        }
        else {
            B = c(1,dbeta(x[k[i]],1:k[i], n0+1-(1:k[i])))
        }
        
        # update step sigma < x[k[i]]
        z = sigma[(n0+1):n]
        for(j in 1:(n-n0)) {
            B[2:(k[i]+1)] = (1-z[j]) * B(2:(k[i]+1)) + z[j] * B[1:k[i]]
        }
        
        Beta[i] = B[k[i]+1]
        
    }
    
    names(Beta) = k
    
    return(Beta)

}

"correctBetaPvalues" = function( p, k ) {
    
    pval = dbeta(p,1,k)
    return(pval)
    
}
