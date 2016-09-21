# perform T-SNE
"tsne" <- function(X, initial_config = NULL, k = 2, max_iter = 1000, min_cost = 0, epoch = 100) {
    
    cat("Performing t-SNE.\n")
    
    momentum = 0.8
    final_momentum = 0.8
    mom_switch_iter = 250

    epsilon = 500
    min_gain = 0.01
    initial_P_gain = 4
    
    n = nrow(X)

    eps = .Machine$double.eps

    if (!is.null(initial_config) && is.matrix(initial_config)) {         
        if (nrow(initial_config) != n | ncol(initial_config) != k) {
            stop('initial_config argument does not match necessary configuration for X')
        }
        ydata = initial_config
        initial_P_gain = 1
        
    }
    else {
        ydata = matrix(rnorm(k * n),n)
    }
    
    P = X
    P = 0.5 * (P + t(P))

    P[P < eps]<-eps
    P = P/sum(P)
    
    P = P * initial_P_gain
    grads = matrix(0,nrow(ydata),ncol(ydata))
    incs = matrix(0,nrow(ydata),ncol(ydata))
    gains = matrix(1,nrow(ydata),ncol(ydata))
    
    for (iter in 1:max_iter) {
        
        if (iter %% epoch == 0) {
            cost = sum(apply(P * log((P+eps)/(Q+eps)),1,sum))
            
            cat("Epoch: Iteration #",iter," error is: ",cost,"\n")
            
            if (cost < min_cost) {
                break
            }
        }
        
        sum_ydata = apply((ydata^2),1,sum)
        num =  1/(1 + sum_ydata + sweep(-2*ydata %*% t(ydata),2,-t(sum_ydata)))
        diag(num) = 0
        Q = num / sum(num)
        
        if (any(is.nan(num))) {
            message ('NaN in grad. descent')
        }
        
        Q[Q < eps] = eps
        
        stiffnesses = (P-Q) * num
        grads = 4 * (diag(apply(stiffnesses,2,sum)) - stiffnesses) %*% ydata
        
        gains = (gains + .2) * abs(sign(grads) != sign(incs)) + gains * .8 * abs(sign(grads) == sign(incs))
        gains[gains < min_gain] = min_gain

        incs = momentum * incs - epsilon * (gains * grads)
        ydata = ydata + incs
        ydata = sweep(ydata,2,apply(ydata,2,mean))
        
        # we are constraining the ydata
        ydata[ydata < -100] = -100
        ydata[ydata > 100] = 100
        
        if (iter == mom_switch_iter) {
            momentum = final_momentum
        }
        
        if (iter == 100 && is.null(initial_config)) {
            P = P/4
        }
        
    }
    
    return(ydata)
    
}
