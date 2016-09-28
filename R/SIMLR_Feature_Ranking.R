# A is the similarity matrix by SIMLR
# X is the data of size nxp
"SIMLR_Feature_Ranking" <- function( A, X ) {
    
    # traspose X
    X = t(X)
    
    # start the computation
    res = lapply(1:100,FUN=function( x ) {
        cat(x,"\n")
        index = sample(1:nrow(A))
        index = index[1:round(nrow(A)*0.9)]
        Ai = A[index,index]
        Xi = X[index,]
        res = LaplacianScore(Xi,Ai)
        return(res)
    })
    yscore = Reduce("rbind",res)

    yscore = 1 - yscore
    numerator = (t(yscore) - min(as.vector(yscore)) + .Machine$double.eps)
    denominator = (max(as.vector(yscore)) - min(as.vector(yscore)) + .Machine$double.eps)
    glist = numerator / denominator
    res = aggregateRanks(glist)
    res2 = sort(res$pval,index.return=TRUE)
    
    res = list(pval=res2$x,aggR=res2$ix)
    
    return(res)
    
}

"LaplacianScore" = function( X, W ) {
    
    nSmp = nrow(X)
    nFea = ncol(X)
    
    if(nrow(W) != nSmp) {
        stop('W is error')
    }

    D = apply(W,MARGIN=1,FUN=sum)
    D = matrix(D,nrow=length(D),ncol=1)
    L = W

    allone = rep(1,nSmp)

    tmp1 = t(D) %*% X
    
    D_temp = array(0,c(length(D),length(D)))
    diag(D_temp) = D
    D = D_temp
    
    DPrime = colSums(t(t(X)%*%D)*X)-tmp1*tmp1/sum(diag(D))
    LPrime = colSums(t((t(X)%*%L))*X)-tmp1*tmp1/sum(diag(D))
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
    
    if(method=="RRA") {
        aggR = rhoScores(rmat,topCutoff)
        pval = aggR
    }
    else {
        stop('Method should be one of:  "min", "geom.mean", "mean", "median" or "RRA"')
    }
    
    return(list(pval=pval,aggR=aggR))
    
}

"rhoScores" = function( r, topCutoff = NaN ) {
    
    rho = rep(NaN,nrow(r))
    
    for(rInd in 1:nrow(r)) {
        r1 = r[rInd,]
        if(is.nan(topCutoff)) {
            x = betaScores(r1)
            rho[rInd] = correctBetaPvalues(min(x),length(which(!is.nan(x))))
        }
    }
    
    return(rho)
    
}

"betaScores" = function( r ) {
    
    n = length(which(!is.nan(r)))
    p = rep(NA,length(r))
    r = sort(r)
    p[1:n] = pbeta(r[1:n],1:n,n:1)
    return(p)
    
}

"correctBetaPvalues" = function( p, k ) {
    
    pval = pbeta(p,1,k)
    return(pval)
    
}
