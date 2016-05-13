
function K=calckernel(kernel_type,kernel_param,X1,X2);
% CALCKERNEL Computes Gram matrix of a specified kernel
% Usage:
% K=calckernel(kernel_type,kernel_param,X1);
% K=calckernel(kernel_type,kernel_param,X1,X2);
% 
% kernel_type = 'linear' | 'poly' | 'rbf'
% kernel_param = -- | degree | sigma
% 
% Given a single data matrix X (n x d where d is dimensionality)
% returns Gram matrix K (n x n)
%
% Given two data matrices X1 (n1 x d), X2 (n2 x d)
% returns Gram matrix K (n1 x n2)
 


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
    
end    
