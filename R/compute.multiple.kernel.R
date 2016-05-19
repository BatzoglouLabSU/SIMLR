# compute and returns the multiple kernel
"multiple.kernel" = function( x ) {
	
	###
	kerneltype={ 'poly' };
	kernelpara={ [0]};
	Kernels=calcmutikernel(kerneltype,kernelpara,x,x);
	N = size(Kernels,1);
	KK = size(Kernels,3);
	sigma = [2:-0.25:1];
	Diff = (dist2(x)).^2;
	[T,INDEX]=sort(Diff,2);
	[m,n]=size(Diff);
	allk = 10:2:30;
	t=1;
	for l = 1:length(allk)
	    l
	    if allk(l) < (size(x,1)-1)
	        TT=mean(T(:,2:(allk(l)+1)),2)+eps;
	        Sig=(repmat(TT,1,n)+repmat(TT',n,1))/2;
	        Sig=Sig.*(Sig>eps)+eps;
	        for j = 1:length(sigma)
	            W=normpdf(Diff,0,sigma(j)*Sig);
	            Kernels(:,:,KK+t) = (W + W')/2;
	            t = t+1;
	        end
	    end
	end
	clear T
	clear INDEX
	clear W
	clear Sig
	for i = 1:size(Kernels,3)
	    K = Kernels(:,:,i);
	    k = 1./sqrt(diag(K)+1);
	    G = K.*(k*k');
	    clear K
	    D_Kernels(:,:,i) = (repmat(diag(G),1,length(G)) +repmat(diag(G)',length(G),1) - 2*G)/2;
	    D_Kernels(:,:,i) = D_Kernels(:,:,i) - diag(diag(D_Kernels(:,:,i)));
	end
	###
	
}

# compute the multiple kernel
"compute.multiple.kernel" = function( kernel.type, kernel.params, x1, x2 = NA ) {
	
	# set the parameters for x1
	n1 = nrow(x1)
	d1 = ncol(x1)
	
	# set the parameters for x2
	if(isna(x2)) {
		n2 = n1
		d2 = d1
		flag = 0
	}
	else {
		n2 = nrow(x2)
		d2 = ncol(x2)
		if(d1!=d2) {
			stop("Error in the data.")
		}
		flag = 1
	}
	
	# set the count.kernel.type
	count.kernel.type = length(kernel.params)
	if(length(kernel.type)<1 || length(kernel.type)!=count.kernel.type) {
		stop("Error in the parameters.")
	}
	
	# count the parameters
	k.count = 0
	for (i in 1:count.kernel.type) {
		k.count = k.count + length(kernel.params[[i]])
	}
	
	# set the structure to save the kernels
	kernels = array(0,c(n1,n2,k.count))
	
	# compute the kernels
	iteration.kernels = 1
	for (i in 1:count.kernel.type) {
		for (j in 1:length(kernel.params[[i]]) {
			single.kernel.type = kernel.type[[i]]
			single.kernel.parameters = kernel.params[[i]][j]
			if(flag==3) {
				kernels(,,iteration.kernels) = compute.kernel(single.kernel.type, single.kernel.parameters,x1)
			}
			else {
				kernels(,,iteration.kernels) = compute.kernel(single.kernel.type, single.kernel.parameters,x1,x2)
			}
		}
		iteration.kernels = iteration.kernels + 1
	}
	
}

# compute the single kernel
"compute.kernel" = function( kernel.type, kernel.params, x1, x2 = NA ) {
	
	###
	dim2=size(X1,2);
	n2=size(X1,1);
	
	
	if nargin==4
	n1=size(X2,1);
	dim=size(X2,2);
	if(dim~=dim2)
	s=sprintf('dim1: %d ,dim2: %d ',dim2,dim);
	error(s);
	end;
	end
	
	switch kernel_type
	
	case 'linear'
	    
	    if nargin==4
	        K=X1*X2';  
	    else
	        K=X2*X2';
	    end
	    
	case 'poly'
	    
	    if nargin==4
	        K=(X1*X2').^kernel_param;
	    else
	        K=(X1*X1').^kernel_param;
	    end
	    
	case 'rbf'  
    
    if nargin==4
        K = exp(-(repmat(sum(X2.*X2,2)',n2,1) + repmat(sum(X1.*X1,2),1,n1) ...
            - 2*X1*X2')/(2*kernel_param^2)); 
    else
    
        P=sum(X1.*X1,2);
        K = exp(-(repmat(P',n2,1) + repmat(P,1,n2) ...
            - 2*X1*X1')/(2*kernel_param^2)); 
        
    end
	###
	
}
