function [y, S, F, ydata,alphaK,timeOurs,converge,LF] = SIMLR_IMPROVED(X, c, no_dim, k, ifimpute,normalize,ytrue)

%%%
tol = 0.0001;
if nargin==2
    no_dim=c;
    k=min(20,round(size(X,2)/c/2));
    ifimpute = 0;
    normalize = 1;
    ytrue=[];
end

if nargin==3
    k=min(20,round(size(X,2)/c/2));
    ifimpute = 0;
    normalize = 1;
    ytrue=[];
end

if nargin==4
    ifimpute = 0;
    normalize = 1;
    ytrue=[];
end


if ifimpute
    X = X';
    [I,J] = find(X==0);
    Xmean = mean(X);
    X(sub2ind(size(X),I,J)) = Xmean(J);
    X = X';
end

normalize=0;
if normalize
    X = X';
    %     X = X - min(X(:));
    %     X = X / max(X(:));
    %     X = bsxfun(@minus, X, mean(X, 1));
    X = normmean0std1(X);
    X = X';
end

t0 = tic;


MAX_NITER = 30;
num = size(X,2);
r = -1;
beta = 0.8;
D_Kernels = multipleK(X');
alphaK = 1/(size(D_Kernels,3))*ones(1,size(D_Kernels,3));
distX = mean(D_Kernels,3);
[distX1, idx] = sort(distX,2);
A = zeros(num);
di = distX1(:,2:(k+2));
rr = 0.5*(k*di(:,k+1)-sum(di(:,1:k),2));
id = idx(:,2:(k+2));
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
S = (1-beta)*S0+beta*A0;
%S= S0;
D0 = diag(sum(S));
L0= D0-S;
[F, temp, evs]=myeig(S, c);
%[F, temp, evs]=eig_fast_bo(X', S, c);
for iter = 1:MAX_NITER
    iter
    distf = L2_distance_1(F',F');
    A = zeros(num);
    %b = idx(:,2:(2*k+2));
    b = idx(:,2:end);
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
    [F, temp, ev]=myeig( S, c);
    evs(:,iter+1) = ev;
    %Dall = (D_Kernels).*repmat(S,1,1,size(D_Kernels,3));
    for i = 1:size(D_Kernels,3)
        temp = (eps+D_Kernels(:,:,i)).*(eps+S);
        DD(i) = mean(sum(temp));
    end
    alphaK0 = umkl_bo(DD);
    alphaK0 = alphaK0/sum(alphaK0);
    alphaK = (1-beta)*alphaK + beta*alphaK0;
    alphaK = alphaK/sum(alphaK);
    fn1 = sum(ev(1:c));
    fn2 = sum(ev(1:c+1));
    converge(iter) = ev(c+1)-ev(c);
    fn2-fn1
    if iter<10
        if (ev(end) > 0.0001)
            lambda = (1.1+0.1*iter)*lambda;
            %r = r/1.1;
        end
    else
        if (converge(iter)>1.01*converge(iter-1))
            S = S_old;
            if converge(iter-1) < 0.001
                warning('Maybe you should set a larger value of c');
            end
            break;
        elseif (abs(converge(iter) - converge(iter-1))<tol)&((abs(converge(iter-1) - converge(iter-2))<tol))
            break;
        else
            lambda = 2*lambda;
        end
    end
    S_old = S;
    distX = Kbeta(D_Kernels,alphaK');
    [distX1, idx] = sort(distX,2);
end;
S = Network_Diffusion(S,2*k);
D = diag(sum(S));
L = D - S;
[F, temp, ev]=myeig( S, c);
LF = F;
y = litekmeans(F, c,'replicates',200);
% if length(no_dim)==1
%     F = tsne_p((S),[], U(:,1:no_dim));
% else
%     clear F;
%     for i = 1:length(no_dim)
%         F{i} = tsne_p((S),[], U(:,1:no_dim(i)));
%     end
% end
timeOurs = toc(t0);

%y = litekmeans(F, c,'replicates',20);
%y=[];

F = tsne_p(S,[],c);
ydata = tsne_p(S);
%ydata = [];
%F=[];
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

function D_Kernels = multipleK(x)

%kerneltype={'poly'};
%kernelpara={[0 0.5 1 2]};
%Kernels=calcmutikernel(kerneltype,kernelpara,x,x);
Kernels = ones(size(x,1),size(x,1),1);
KK = size(Kernels,3);
Diff =sqrt(dist2(x));
[T,INDEX]=sort(Diff,2);
[m,n]=size(Diff);
sigma = [2:0.25:2.5];
allk = 10:2:20;
%sigma = [2];
%allk = 10:10:20;
t=1;
for l = 1:length(allk)
    l
    if allk(l) < (size(x,1)-1)
        TT=mean(T(:,2:(allk(l)+1)),2)+eps;
        Sig=(repmat(TT,1,n)+repmat(TT',n,1))/2;
        Sig=Sig.*(Sig>eps)+eps;
        for j = 1:length(sigma)
            W=normpdf_bo(Diff,0,sigma(j)*Sig);
            Kernels(:,:,KK+t) = (W + W')/2;
            t = t+1;
        end
    end
end

for i = 1:size(Kernels,3)
    K = Kernels(:,:,i);
    k = 1./sqrt(sum(K.^2)+1);
    G = K.*(k*k');
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
W = (dominateset(double(abs(A)),min(2*K,length(A)-1))).*sign(A);
A = A-diag(diag(A));
P = (dominateset(double(abs(A)),min(K,length(A)-1))).*sign(A);
DD = sum(abs(P'));
P = P + (eye(length(P))+diag(sum(abs(P'))));
P = dn(P,'ave');
% [U,D] = eigs(sparse(P),K,'LM');
% d = real((diag(D))+eps);
% d = d./max(d);
% alpha = 0.8;
% beta = 1;
% d = ((1-alpha)*d.^beta./(1-alpha*d.^beta));
% D = diag(real(d));
% W = U*D*U';
% W = W + P;

W = sparse(W);P = sparse(P);
for i = 1:3
   % W = mtimesx_sparse(W,'N',P,'N') + eye(length(W));
W = W*P + eye(length(W));
end
W = W+W';
%W = (dominateset(double(abs(W)),min(3*K,length(A)-1))).*sign(A);
W = dn(W,'ave');

W = (W.*(1-eye(length(W))))./repmat(1-diag(W),1,length(W));
D=sparse(1:length(DD),1:length(DD),DD);
W=D*(W);
W = (W+W')/2;
end

function y = normpdf_bo(x,mu,sigma)

y = exp(-0.5 * ((x - mu)./sigma).^2)./ (sqrt(2*pi) .* sigma);
end

function [eigvec, eigval, eigval_full] = myeig(S, c)

D0 = diag(sum(S));
A= D0-S;
%[v d] = eig(A);
%d = diag(d);
OPTS.issym = 1;
%[v,d] = laneig(sparse(real(A+A'+eps)/2),c+1,'SM',OPTS);
%[v,d] = eigs(sparse(real(A)),c+1,'SM',OPTS);
%tic;[v,d] = laneig((real(A+A'+eps)/2),c+1,'SM',OPTS);toc
tic;[v,d] = eigs(sparse(real(A)),c+1,'SM',OPTS);toc
d =  diag(real(d));



%d = real(d);
[d1, idx] = sort(d,'ascend');
idx1 = idx(1:c);
eigval = d(idx1);
%eigvec = real(v(:,idx1)*diag(sqrt(1-eigval)));
eigvec = real(v(:,idx1));
eigval_full = d(idx);
end












