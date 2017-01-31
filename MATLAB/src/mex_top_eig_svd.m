function [U,d] = mex_top_eig_svd(S0, ind, K)

I = repmat([1:size(S0,1)]',1,size(S0,2));
S0 = sparse(I(:), ind(:), S0(:), size(S0,1),size(S0,1)) + sparse(ind(:), I(:), S0(:), size(S0,1),size(S0,1));
[U, d, ~] = eigs(S0, K);
d= sqrt(d);
%U = U(:,1:K)*diag(sqrt(sum(U(:,1:K).^2)));
