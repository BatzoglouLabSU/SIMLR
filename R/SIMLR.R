# perform the SIMLR clustering algorithm
"SIMLR" <- function( X, c , no.dim = NA, k = 10, if.impute = FALSE, normalize = FALSE ) {
    
    # set any required parameter to the defaults
    if(is.na(no.dim)) {
        no.dim = c
    }
    
    # check the if.impute parameter
    if(if.impute == TRUE) {
        X = t(X)
        X_zeros = which(X==0,arr.ind=TRUE)
        if(length(X_zeros)>0) {
            R_zeros = as.vector(X_zeros[,"row"])
            C_zeros = as.vector(X_zeros[,"col"])
            ind = (C_zeros - 1) * nrow(X) + R_zeros
            X[ind] = as.vector(colMeans(X))[C_zeros]
        }
        X = t(X)
    }
    
    # check the normalize parameter
    if(normalize == TRUE) {
        X = t(X)
        X = X - min(as.vector(X))
        X = X / max(as.vector(X))
        C_mean = as.vector(colMeans(X))
        X = apply(X,MARGIN=1,FUN=function(x) return(x-C_mean))
    }
    
    # start the clock to measure the execution time
    ptm = proc.time()
    
    # set some parameters
    NITER = 30
    num = ncol(X)
    r = -1
    beta = 0.8
    
    cat("Computing the Kernels...")
    
    # compute the kernels
    D_Kernels = multiple.kernel(t(X))
    
    # set up some parameters
    alphaK = 1 / rep(dim(D_Kernels)[3],dim(D_Kernels)[3])
    distX = t(apply(D_Kernels,MARGIN=3,FUN=mean))
    
    # sort distX for rows
    res = apply(distX,MARGIN=1,FUN=function(x) return(sort(x,index.return = TRUE)))
    distX1 = array(0,c(nrow(distX),ncol(distX)))
    idx = array(0,c(nrow(distX),ncol(distX)))
    for(i in 1:nrow(distX)) {
        distX1[i,] = res[[i]]$x
        idx[i,] = res[[i]]$ix
    }
    
    A = rep(0,num)
    di = distX1[,2:(k+2)]
    rr = 0.5 * (k*di[k+1] - sum(di[1:k]))
    id = idx[,2:k+2]
    
    numerator = (apply(array(0,c(length(di),length(di))),MARGIN=2,FUN=function(x) {x=di}) - di)
    temp = (k*di[k+1] - sum(di[k+1]) + .Machine$double.eps)
    denominator = apply(array(0,c(length(tmp),length(tmp))),MARGIN=1,FUN=function(x) {x=tmp})
    temp = numerator / denominator
    a = apply(array(0,c(length(t(1:num)),length(t(1:num)))),MARGIN=2,FUN=function(x) {x=1:num})
    A[(id-1)*length(A) + a] = temp
    
    if(r<=0) {
        r = mean(rr)
    }
    lambda = max(mean(rr),0)
    A[is.nan(A)] = 0
    A0 = (A + t(A)) / 2
    S0 = 1 - distX
    
    cat("Performing network diffiusion...")
    
    # perform network diffiusion
    S0 = network.diffusion(S0,k)
    
    # compute dn
    S0 = dn(S0,'gph')
    S = (1 - beta) %*% S0 + beta %*% A0
    D0 = diag(sum(S))
    L0 = D0 - S
    
    # ####### TO DO
    
    # [F, temp, evs]=eig1(L0, c, 0)
    
    # for iter = 1:NITER
        # distf = L2_distance_1(F',F');
        # A = zeros(num);
        # b = idx(:,2:(2*k+2));
        # a = repmat([1:num]',1,size(b,2));
        # inda = sub2ind(size(A),a(:),b(:));
        # ad = reshape((distX(inda)+lambda*distf(inda))/2/r,num,size(b,2));
        # ad = projsplx_c(-ad')';
        # A(inda) = ad(:);
        # A(isnan(A))=0;
        # A = (A+A')/2;
        # S = (1-beta)*S+beta*A;
        # S = Network_Diffusion(S,k);
        # D = diag(sum(S));
        # L = D - S;
        # F_old = F;
        # [F, temp, ev]=eig1(L, c, 0);
        # evs(:,iter+1) = ev;
        # for i = 1:size(D_Kernels,3)
        # temp = D_Kernels(:,:,i).*S;
            # DD(i) = mean(sum(temp-diag(diag(temp))));
        # end
        # alphaK0 = umkl_bo(DD);
        # alphaK0 = alphaK0/sum(alphaK0);
        # alphaK = (1-beta)*alphaK + beta*alphaK0;
        # alphaK = alphaK/sum(alphaK);
        # fn1 = sum(ev(1:c));
        # fn2 = sum(ev(1:c+1));
        # converge(iter) = fn2-fn1;
        # fn2-fn1
        # if iter<10
            # if (ev(end) > 0.000001)
                # lambda = 1.5*lambda;
                # r = r/1.01;
            # end
        # else
            # if (converge(iter)>converge(iter-1))
                # S = S_old;
                # if converge(iter-1) > 0.2
                    # warning('Maybe you should set a larger value of c');
                # end
                # break;
            # end
        # end
        # S_old = S;
        # distX = Kbeta(D_Kernels,alphaK');
        # [distX1, idx] = sort(distX,2);
    # end;
    # LF = F;
    # S = Network_Diffusion(S,2*k);
    
    # D = diag(sum(S));
    # L = D - S;
    # [U,D] = eig(L);
    # if length(no_dim)==1
        # F = tsne_p((S),[], U(:,1:no_dim));
    # else
        # clear F;
        # for i = 1:length(no_dim)
            # F{i} = tsne_p((S),[], U(:,1:no_dim(i)));
        # end
    # end
    
  # '  
      
    # ####### TO DO
    
    # compute the execution time
    execution.time = proc.time() - ptm
    
}
