% Kbeta
%
% Computes weighted sum of kernel matrices
%
% Usage: 
%	K = Kbeta(Ks,w);
%	K = Kbeta(Ks,w,symmetric);
%
% Input:
%	Ks - Gram matrices [n x n x M]
%	w - weights [M x 1]
%	symmetric (=false) - optional argument. If 1 the Ks
%		matrices are assumed to be symmetric which results in a
%		small speedup. 
%
% Output: 
%	K - weighted kernel matrix equiv with
%		K=0;for m=1:M, K=K+ w(m)*Ks(:,:,m);end
%
% Example: 
%	
%	K = Kbeta(randn(100,100,10),rand(10,1));
%
% Peter Gehler 07/2008 pgehler@tuebingen.mpg.de
