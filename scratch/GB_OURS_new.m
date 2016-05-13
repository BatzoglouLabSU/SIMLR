function [y, S, F, ydata,alphaK,timeOurs,converge,LF] = GB_OURS_new(X, c,no_dim, k, ifimpute,normalize,ytrue)

%%%
if nargin==2
    no_dim=c;
    k=20;
    ifimpute = 0;
    normalize = 0;
    ytrue=[];
end

if nargin==3
    k=20;
    ifimpute = 0;
    normalize = 0;
    ytrue=[];
end

if nargin==4
    ifimpute = 0;
    normalize = 0;
    ytrue=[];
end

if ifimpute
    X = X';
    [I,J] = find(X==0);
    Xmean = mean(X);
    X(sub2ind(size(X),I,J)) = Xmean(J);
    X = X';
end

if normalize
    X = X';
    X = X - min(X(:));
    X = X / max(X(:));
    X = bsxfun(@minus, X, mean(X, 1));
    X = X';
end

t0 = tic;


NITER = 30;
num = size(X,2);
r = -1;
beta = 0.8;
D_Kernels = multipleK(X',c);
alphaK = 1/(size(D_Kernels,3))*ones(1,size(D_Kernels,3));
distX = mean(D_Kernels,3);
[distX1, idx] = sort(distX,2);
A = zeros(num);
di = distX1(:,2:(k+2));
rr = 0.5*(k*di(:,k+1)-sum(di(:,1:k),2));
id = idx(:,2:k+2);
temp = (repmat(di(:,k+1),1,size(di,2))-di)./repmat((k*di(:,k+1)-sum(di(:,1:k),2)+eps),1,size(di,2));
a = repmat([1:num]',1,size(id,2));
A(sub2ind(size(A),a(:),id(:)))=temp(:);
if r <= 0
    r = mean(rr);
end
lambda = max((mean(rr)),0);
A(isnan(A))=0;
A0 = (A+A')/2;
S0 = 1-distX;
S0 = Network_Diffusion(S0,k);
S0 = dn(S0,'gph');
%S = (1-beta)*S0+beta*A0;
%S= S0;
S = (1-beta)*S0+beta*A0;
D0 = diag(sum(S));
L0= D0-S;
[F, temp, evs]=eig1(L0, c, 0);

for iter = 1:NITER
    distf = L2_distance_1(F',F');
    A = zeros(num);
    b = idx(:,2:(2*k+2));
    a = repmat([1:num]',1,size(b,2));
    inda = sub2ind(size(A),a(:),b(:));
    ad = reshape((distX(inda)+lambda*distf(inda))/2/r,num,size(b,2));
    ad = projsplx_c(-ad')';
    A(inda) = ad(:);
    A(isnan(A))=0;
    A = (A+A')/2;
    S = (1-beta)*S+beta*A;
    S = Network_Diffusion(S,k);
    D = diag(sum(S));
    L = D - S;
    F_old = F;
    [F, temp, ev]=eig1(L, c, 0);
    evs(:,iter+1) = ev;
%    Dall = (D_Kernels).*repmat(S,1,1,size(D_Kernels,3));
    for i = 1:size(D_Kernels,3)
	temp = D_Kernels(:,:,i).*S;
        DD(i) = mean(sum(temp-diag(diag(temp))));
    end
    alphaK0 = umkl_bo(DD);
    alphaK0 = alphaK0/sum(alphaK0);
    alphaK = (1-beta)*alphaK + beta*alphaK0;
    alphaK = alphaK/sum(alphaK);
    fn1 = sum(ev(1:c));
    fn2 = sum(ev(1:c+1));
    converge(iter) = fn2-fn1;
    fn2-fn1
    if iter<10
        if (ev(end) > 0.000001)
            lambda = 1.5*lambda;
            r = r/1.1;
        end
    else
        if (converge(iter)>converge(iter-1))
            S = S_old;
            if converge(iter-1) > 0.2
                warning('Maybe you should set a larger value of c');
            end
            break;
        end
    end
    S_old = S;
    distX = Kbeta(D_Kernels,alphaK');
    [distX1, idx] = sort(distX,2);
end;
LF = F;
S = Network_Diffusion(S,2*k);

D = diag(sum(S));
L = D - S;
[U,D] = eig(L);
if length(no_dim)==1
    F = tsne_p((S),[], U(:,1:no_dim));
else
    clear F;
    for i = 1:length(no_dim)
        F{i} = tsne_p((S),[], U(:,1:no_dim(i)));
    end
end
timeOurs = toc(t0);

y = litekmeans(F, c,'replicates',200);
%y=[];

ydata = tsne_p(S);

end
function thisP = umkl_bo(D,beta)
if nargin<2
    beta = 1/length(D);
end
tol = 1e-4;
u = 20;logU = log(u);
[H, thisP] = Hbeta(D, beta);
betamin = -Inf;
betamax = Inf;
% Evaluate whether the perplexity is within tolerance
Hdiff = H - logU;
tries = 0;
while (abs(Hdiff) > tol) && (tries < 30)
    
    % If not, increase or decrease precision
    if Hdiff > 0
        betamin = beta;
        if isinf(betamax)
            beta = beta * 2;
        else
            beta = (beta + betamax) / 2;
        end
    else
        betamax = beta;
        if isinf(betamin)
            beta = beta / 2;
        else
            beta = (beta + betamin) / 2;
        end
    end
    
    % Recompute the values
    [H, thisP] = Hbeta(D, beta);
    Hdiff = H - logU;
    tries = tries + 1;
end
end

function [H, P] = Hbeta(D, beta)
D = (D-min(D))/(max(D) - min(D)+eps);
P = exp(-D * beta);
sumP = sum(P);
H = log(sumP) + beta * sum(D .* P) / sumP;
P = P / sumP;
end

function D_Kernels = multipleK(x,c)

kerneltype={ 'poly' };
kernelpara={ [0]};
Kernels=calcmutikernel(kerneltype,kernelpara,x,x);
N = size(Kernels,1);
KK = size(Kernels,3);
sigma = [2:-0.5:1];
Diff = (dist2(x));
[T,INDEX]=sort(Diff,2);
[m,n]=size(Diff);
allk = 10:5:25;
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

KK = size(Kernels,3);

Diff = (pdist2(x,x,'cosine'));
[T,INDEX]=sort(Diff,2);
[m,n]=size(Diff);
allk = 10:5:25;
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

x0 = x;
%%%%
KK = size(Kernels,3);
x = compute_mapping(x0,'PCA', c + 1);
Diff = (dist2(x));
[T,INDEX]=sort(Diff,2);
[m,n]=size(Diff);
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
KK = size(Kernels,3);
Diff = (pdist2(x,x,'cosine'));
[T,INDEX]=sort(Diff,2);
[m,n]=size(Diff);
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


KK = size(Kernels,3);
x = compute_mapping(x0,'MDS', c + 1);
Diff = (dist2(x));
[T,INDEX]=sort(Diff,2);
[m,n]=size(Diff);
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
KK = size(Kernels,3);
Diff = (pdist2(x,x,'cosine'));
[T,INDEX]=sort(Diff,2);
[m,n]=size(Diff);
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
    %G = K;
    %G = G - diag(diag(G));
    %G = dn(G+1*eye(length(G))*eps,'gph');
    %G = dn(G,'gph');
    D_Kernels(:,:,i) = (repmat(diag(G),1,length(G)) +repmat(diag(G)',length(G),1) - 2*G)/2;
    D_Kernels(:,:,i) = D_Kernels(:,:,i) - diag(diag(D_Kernels(:,:,i)));
    %G = dn(G+eye(length(G))*eps,'gph');
    %D_Kernels(:,:,2*i) = 2 - 2*G - 2*eye(length(G));
end
D_Kernels(:,:,1)=[];
end



%
function W = Network_Diffusion(A, K)
A = A-diag(diag(A));
P = (dominateset(double(abs(A)),min(K,length(A)-1))).*sign(A);
DD = sum(abs(P'));
P = P + (eye(length(P))+diag(sum(abs(P'))));
P = (TransitionFields(P));
[U,D] = eig(P);
d = real((diag(D))+eps);
alpha = 0.5;
beta = 1;
d = (1-alpha)*d.^beta./(1-alpha*d.^beta);


D = diag(real(d));
W = U*D*U';

W = (W.*(1-eye(length(W))))./repmat(1-diag(W),1,length(W));
D=sparse(1:length(DD),1:length(DD),DD);
W=D*(W);
W = (W+W')/2;

end
