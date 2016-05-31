.Hbeta <-
function(D, beta){
	P = exp(-D * beta)
	sumP = sum(P)
	if (sumP == 0){
		H = 0
		P = D * 0
	} else {
		H = log(sumP) + beta * sum(D %*% P) /sumP
		P = P/sumP
	}
	r = {}
	r$H = H
	r$P = P
	r
}

.x2p <-
function(X,perplexity = 15,tol = 1e-5){
	if (class(X) == 'dist') {
		D = X
		n = attr(D,'Size')
	} else{
		D = dist(X)
		n = attr(D,'Size')
	}

	D = as.matrix(D)
	P = matrix(0, n, n )		
	beta = rep(1, n)
	logU = log(perplexity)
	
	for (i in 1:n){
		betamin = -Inf
		betamax = Inf
		Di = D[i, -i]
		hbeta = .Hbeta(Di, beta[i])
		H = hbeta$H; 
		thisP = hbeta$P
		Hdiff = H - logU;
		tries = 0;

		while(abs(Hdiff) > tol && tries < 50){
			if (Hdiff > 0){
				betamin = beta[i]
				if (is.infinite(betamax)) beta[i] = beta[i] * 2
				else beta[i] = (beta[i] + betamax)/2
			} else{
				betamax = beta[i]
				if (is.infinite(betamin))  beta[i] = beta[i]/ 2
				else beta[i] = ( beta[i] + betamin) / 2
			}
			
			hbeta = .Hbeta(Di, beta[i])
			H = hbeta$H
			thisP = hbeta$P
			Hdiff = H - logU
			tries = tries + 1
		}	
			P[i,-i]  = thisP	
	}	
	
	r = {}
	r$P = P
	r$beta = beta
	sigma = sqrt(1/beta)
	
	message('sigma summary: ', paste(names(summary(sigma)),':',summary(sigma),'|',collapse=''))

	r 
}

.whiten <-
function(X, row.norm=FALSE, verbose=FALSE, n.comp=ncol(X))
{  
	n.comp; # forces an eval/save of n.comp
	if (verbose) message("Centering")
   n = nrow(X)
	p = ncol(X)
	X <- scale(X, scale = FALSE)
   X <- if (row.norm) 
       t(scale(X, scale = row.norm))
   else t(X)

   if (verbose) message("Whitening")
   V <- X %*% t(X)/n
   s <- La.svd(V)
   D <- diag(c(1/sqrt(s$d)))
   K <- D %*% t(s$u)
   K <- matrix(K[1:n.comp, ], n.comp, p)
   X = t(K %*% X)
	X
}

