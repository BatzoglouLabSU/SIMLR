# perform network diffusion of K steps over the network A
"network.diffusion" <- function( A, K ) {
    
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
    W = diag(DD) %*% W
    W = (W + t(W)) / 2
    
    W[which(W<0,arr.ind=TRUE)] = 0
    
    return(W)
    
}

# compute the dominate set for the matrix aff.matrix and NR.OF.KNN
"dominate.set" <- function( aff.matrix, NR.OF.KNN ) {
    
    # create the structure to save the results
    PNN.matrix = array(0,c(nrow(aff.matrix),ncol(aff.matrix)))
    
    # sort each row of aff.matrix in descending order and saves the sorted 
    # array and a collection of vectors with the original indices
    res.sort = apply(t(aff.matrix),MARGIN=2,FUN=function(x) {return(sort(x, decreasing = TRUE, index.return = TRUE))})
    sorted.aff.matrix = t(apply(as.matrix(1:length(res.sort)),MARGIN=1,function(x) { return(res.sort[[x]]$x) }))
    sorted.indices = t(apply(as.matrix(1:length(res.sort)),MARGIN=1,function(x) { return(res.sort[[x]]$ix) }))
    
    # get the first NR.OF.KNN columns of the sorted array
    res = sorted.aff.matrix[,1:NR.OF.KNN]
    
    # create a matrix of NR.OF.KNN columns by binding vectors of 
    # integers from 1 to the number of rows/columns of aff.matrix
    inds = array(0,c(nrow(aff.matrix),NR.OF.KNN))
    inds = apply(inds,MARGIN=2,FUN=function(x) {x=1:nrow(aff.matrix)})
    
    # get the first NR.OF.KNN columns of the indices of aff.matrix
    loc = sorted.indices[,1:NR.OF.KNN]
    
    # assign to PNN.matrix the sorted indices
    PNN.matrix[(as.vector(loc)-1)*nrow(aff.matrix)+as.vector(inds)] = as.vector(res)
    
    # compute the final results and return them
    PNN.matrix = (PNN.matrix + t(PNN.matrix))/2
    
    return(PNN.matrix)
    
}

# compute the transition field of the given matrix
"transition.fields" <- function( W ) {
    
    # get any index of columns with all 0s
    zero.index = which(apply(W,MARGIN=1,FUN=sum)==0)
    
    # compute the transition fields
    W = dn(W,'ave')
    
    w = sqrt(apply(abs(W),MARGIN=2,FUN=sum)+.Machine$double.eps)
    W = W / t(apply(array(0,c(nrow(W),ncol(W))),MARGIN=2,FUN=function(x) {x=w}))
    W = W %*% t(W)
    
    # set to 0 the elements of zero.index
    W[zero.index,] = 0
    W[,zero.index] = 0
    
    return(W)
    
}

# normalizes a symmetric kernel
"dn" = function( w, type ) {
    
    # compute the sum of any column
    D = apply(w,MARGIN=2,FUN=sum)
    
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
