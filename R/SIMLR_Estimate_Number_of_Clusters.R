# Estimates the number of clusters by means of two huristics
"SIMLR_Estimate_Number_of_Clusters" = function( X, NUMC = 2:5, cores.ratio = 1 ) {

    D_Kernels = multiple.kernel.numc(t(X),cores.ratio)
    distX = array(0,c(dim(D_Kernels[[1]])[1],dim(D_Kernels[[1]])[2]))
    for (i in 1:length(D_Kernels)) {
        distX = distX + D_Kernels[[i]]
    }
    distX = distX / length(D_Kernels)
    W =  max(max(distX)) - distX
    W = network.diffusion.numc(W,max(ceiling(ncol(X)/20),10))

    Quality = Estimate_Number_of_Clusters_given_graph(W,NUMC)
    Quality_plus = Estimate_Number_of_Clusters_given_graph(W,NUMC+1)
    Quality_minus = Estimate_Number_of_Clusters_given_graph(W,NUMC-1)

    K1 = 2*(1 + Quality) - (2 + Quality_plus + Quality_minus)
    K2 = K1*(NUMC+1)/(NUMC)

    return(list(K1=K1,K2=K2))

}

# This function estimates the number of clusters given the two huristics 
# given in the supplementary materials of our nature method paper. 
# W is the similarity graph; NUMC is a vector which contains the possible choices 
# of number of clusters 
"Estimate_Number_of_Clusters_given_graph" = function( W, NUMC = 2:5 ) {

    quality = NULL

    if (min(NUMC)<1) {
        stop('Note that we always assume a minimum of at least 2 clusters: please change values for NUMC.')
    }

    W = (W + t(W)) / 2

    if(!is.null(NUMC)) {

        degs = rowSums(W)
        D = Matrix(0,nrow=length(degs),ncol=length(degs),sparse=TRUE)
        diag(D) = degs
        L = D - W
        degs[which(degs==0)] = .Machine$double.eps
        # calculate D^(-1/2)
        diag(D) = 1/(degs^0.5)
        # calculate normalized Laplacian
        L = D %*% L %*% D
        # compute the eigenvectors corresponding to the k smallest  eigenvalues
        res = eigen(L)
        eigenvalue = res$values
        U = res$vectors
        res = sort(eigenvalue,decreasing=FALSE,index.return=TRUE)
        eigenvalue = res$x
        b = res$ix
        U = U[,b]
        for (ck in NUMC) {
            if(ck==1) {
                tmp = array(0,dim(U))
                diag(tmp) = 1/(U[,1]+.Machine$double.eps)
                res = sum(sum(tmp%*%U[,1]))
            }
            else {
                UU = U[,1:ck]
                tmp = sqrt(rowSums(UU^2))+.Machine$double.eps
                tmp = matrix(rep(tmp,ck),nrow=length(tmp),ncol=ck)
                UU = UU / tmp
                res = discretisation(UU)
                EigenVectors = res$EigenVectors^2
                temp1 = t(apply(EigenVectors,1,function(x) return(sort(x,decreasing=TRUE))))
                tmp = 1 / (temp1[,1]+.Machine$double.eps)
                tmp1 = Matrix(0,nrow=length(tmp),ncol=length(tmp),sparse=TRUE)
                diag(tmp1) = tmp
                tmp = tmp1%*%temp1[,1:max(2,ck-1)]
                tmp = sum(sum(tmp))
                res = (1-eigenvalue[ck+1])/(1-eigenvalue[ck])*tmp
            }
            quality = c(quality,res)
        }
    
    }

    return(quality)

}

"discretisation" = function( EigenVectors ) {

    n = nrow(EigenVectors)
    k = ncol(EigenVectors)
    vm = sqrt(rowSums(EigenVectors*EigenVectors,2))
    tmp = matrix(rep(vm+.Machine$double.eps,k),nrow=length(vm+.Machine$double.eps),ncol=k)
    EigenVectors = EigenVectors/tmp
    R = array(0,c(k,k))
    R[,1] = t(EigenVectors[round(n/2),])

    c = array(0,c(n,1))
    for(j in 2:k) {
        c = c + abs(EigenVectors%*%R[,(j-1)])
        R[,j] = t(EigenVectors[which.min(c),])
    }

    lastObjectiveValue = 0
    exitLoop = 0
    nbIterationsDiscretisation = 0
    nbIterationsDiscretisationMax = 20
    while(exitLoop==0) {
        nbIterationsDiscretisation = nbIterationsDiscretisation + 1
        EigenvectorsDiscrete = discretisationEigenVectorData(EigenVectors%*%R)
        res_svd = svd(t(EigenvectorsDiscrete)%*%EigenVectors+.Machine$double.eps)
        U = res_svd$u
        S = res_svd$d
        V = res_svd$v
        NcutValue=2*(n-sum(S))
        if(abs(NcutValue-lastObjectiveValue) < .Machine$double.eps || nbIterationsDiscretisation > nbIterationsDiscretisationMax) {
            exitLoop = 1
        }
        else {
            lastObjectiveValue = NcutValue
            R = V%*%t(U)
        }
    }

    res = list(EigenvectorsDiscrete=EigenvectorsDiscrete,EigenVectors=EigenVectors)
    return(res)

}

# Discretizes previously rotated eigenvectors in discretisation 
# Timothee Cour, Stella Yu, Jianbo Shi, 2004 
"discretisationEigenVectorData" = function( EigenVector ) {

    n = nrow(EigenVector)
    k = ncol(EigenVector)

    J = apply(t(EigenVector),2,function(x) return(which.max(x)))
    Y = sparseMatrix(i=1:n,j=t(J),x=1,dims=c(n,k))

    return(Y)

}

# compute and returns the multiple kernel
"multiple.kernel.numc" = function( x, cores.ratio = 1 ) {
    
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
    Diff = dist2(x)
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
        k = 1/sqrt(diag(K)+.Machine$double.eps)
        G = K * (k %*% t(k))
        G1 = apply(array(0,c(length(diag(G)),length(diag(G)))),MARGIN=2,FUN=function(x) {x=diag(G)})
        G2 = t(G1)
        D_Kernels_tmp = (G1 + G2 - 2*G)/2
        D_Kernels_tmp = D_Kernels_tmp - diag(diag(D_Kernels_tmp))
        D_Kernels[[i]] = Matrix(D_Kernels_tmp, sparse=TRUE, doDiag=FALSE)
    }
    
    return(D_Kernels)
    
}

# perform network diffusion of K steps over the network A
"network.diffusion.numc" <- function( A, K ) {
    
    # set the values of the diagonal of A to 0
    diag(A) = 0
    
    # compute the sign matrix of A
    sign_A = A
    sign_A[which(A>0,arr.ind=TRUE)] = 1
    sign_A[which(A<0,arr.ind=TRUE)] = -1
    
    # compute the dominate set for A and K
    P = dominate.set(abs(A),min(K,nrow(A)-1)) * sign_A
    
    # sum the absolute value of each row of P
    DD = apply(abs(P),MARGIN=1,FUN=sum)
    
    # set DD+1 to the diagonal of P
    diag(P) = DD + 1
    
    # compute the transition field of P
    P = transition.fields(P)
    
    # compute the eigenvalues and eigenvectors of P
    eigen_P = eigen(P)
    U = eigen_P$vectors
    D = eigen_P$values
    
    # set to d the real part of the diagonal of D
    d = Re(D + .Machine$double.eps)
    
    # perform the diffusion
    alpha = 0.8
    beta = 2
    d = ((1-alpha)*d)/(1-alpha*d^beta)
    
    # set to D the real part of the diagonal of d
    D = array(0,c(length(Re(d)),length(Re(d))))
    diag(D) = Re(d)
    
    # finally compute W
    
    W = U %*% D %*% t(U)
    diagonal_matrix = array(0,c(nrow(W),ncol(W)))
    diag(diagonal_matrix) = 1
    W = (W * (1-diagonal_matrix)) / apply(array(0,c(nrow(W),ncol(W))),MARGIN=2,FUN=function(x) {x=(1-diag(W))})
    diag(D) = diag(D)[length(diag(D)):1]
    W = (W + t(W)) / 2
    
    W[which(W<0,arr.ind=TRUE)] = 0
    
    return(W)
    
}
