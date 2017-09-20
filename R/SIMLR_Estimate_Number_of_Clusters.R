# Estimates the number of clusters by means of two huristics
"SIMLR_Estimate_Number_of_Clusters" = function( X, NUMC = 2:5, cores.ratio = 1 ) {

    D_Kernels = multiple.kernel(X,cores.ratio)
    distX = array(0,c(dim(D_Kernels[[1]])[1],dim(D_Kernels[[1]])[2]))
    for (i in 1:length(D_Kernels)) {
        distX = distX + D_Kernels[[i]]
    }
    distX = distX / length(D_Kernels)
    W =  max(max(distX)) - distX
    W = network.diffusion(W,max(ceiling(nrow(X)/20),10))

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
