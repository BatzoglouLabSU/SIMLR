function X = fast_pca(in_X, K)

in_X = in_X - repmat(mean(in_X),size(in_X,1),1);
[U, S, ~] = rsvd(in_X, K);
K = min(size(S,2),K);
X = U(:,1:K)*diag(sqrt(diag(S(1:K,1:K))));
X = X./repmat(sqrt(sum(X.*X,2)),1,K);
end

function [U,S,V] = rsvd(A,K)
%-------------------------------------------------------------------------------------
% random SVD
% Extremely fast computation of the truncated Singular Value Decomposition, using
% randomized algorithms as described in Halko et al. 'finding structure with randomness
%
% usage : 
%
%  input:
%  * A : matrix whose SVD we want
%  * K : number of components to keep
%
%  output:
%  * U,S,V : classical output as the builtin svd matlab function
%-------------------------------------------------------------------------------------
% Antoine Liutkus  (c) Inria 2014

[M,N] = size(A);
P = min(2*K,N);
X = randn(N,P);
Y = A*X;
W1 = orth(Y);
B = W1'*A;
[W2,S,V] = svd(B,'econ');
U = W1*W2;
K=min(K,size(U,2));
U = U(:,1:K);
S = S(1:K,1:K);
V=V(:,1:K);
end