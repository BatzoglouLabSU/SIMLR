function [S0, F] = SIMLR_Large_Scale(X, c, k, ifimpute,normalize)

%%%
if nargin==2
    k=10;
    ifimpute = 0;
    normalize = 0;
end

if nargin==3
    ifimpute = 0;
    normalize = 0;
end

if nargin==4
    normalize = 0;
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



NITER = 5;
r = -1;
beta = 0.8;

[ind, val] = mex_KNN_Annoy(X,2*k);
ind = ind + 1;
val = double((abs(val)));
%%%%construct multiple kernels
D_Kernels = mex_multipleK(val, ind,k); %D_Kernels are of size Nx(2k)x55
clear val
alphaK = 1/(size(D_Kernels,3))*ones(1,size(D_Kernels,3));
distX = mean(D_Kernels,3);

di = distX(:,2:(k+2));
rr = 0.5*(k*di(:,k+1)-sum(di(:,1:k),2));

if r <= 0
    r = mean(rr);
end
lambda = max((mean(rr)),0);
clear rr di
%A(isnan(A))=0;
S0 = max(max(distX))-distX;
%S0(:,1) = S0(:,1) + sum(S0,2);


S0 = NE_dn(S0,'ave');


[F,evalues] = mex_top_eig(S0, ind, c);
d = evalues/max(evalues);
d = (1-beta)*d./(1-beta*d.^2);
F =F.*repmat((d+eps)', length(F),1);

F = NE_dn(F,'ave');
F0 = F;
for iter = 1:NITER
    distf = mex_L2_distance_1(F, ind);
    distf = (distX+lambda*distf)/2/r;
    distf = projsplx_c(-distf')';
    S0 = (1-beta)*S0+(beta)*distf;
    [F,evalues] = mex_top_eig(S0, ind, c);
    d = evalues/max(evalues);
    d = (1-beta)*d./(1-beta*d.^2);
    F =F.*repmat((d+eps)', length(F),1);

    F = NE_dn(F,'ave');
    F = (1-beta)*F0+(beta)*F;
    F0 = F;
    for i = 1:size(D_Kernels,3)
        temp = (eps+D_Kernels(:,:,i)).*(eps+S0);
        DD(i) = mean(sum(temp));
    end
    alphaK0 = umkl_bo(DD);
    alphaK0 = alphaK0/sum(alphaK0);
    alphaK = (1-beta)*alphaK + beta*alphaK0;
    alphaK = alphaK/sum(alphaK);
    lambda = 1.5*lambda;
    r = r/1.1;
    distX = Kbeta(D_Kernels,alphaK');
end;
val = S0;

I = repmat([1:size(S0,1)]',1,size(S0,2));
S0 = sparse(I(:), ind(:), S0(:)/2.0, size(S0,1),size(S0,1)) + sparse(ind(:), I(:), S0(:)/2.0, size(S0,1),size(S0,1));
ydata = [];
y = [];
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



function D_kernels = mex_multipleK(val,ind,KK)
if nargin<3
    KK=20;
end

val = val.*val;
sigma = [2:-0.25:1];
allk = [ceil(KK/2):ceil(KK/10):(ceil(KK*1.5))];
t = 1;
for i = 1:length(allk)
    if allk(i)< size(val,2)
    temp = mean(val(:,1:allk(i)),2);
    temp0 = .5*(repmat(temp,1,size(val,2)) + temp(ind))+eps;
    
    for j = 1:length(sigma)
        temp = normpdf(val, 0, sigma(j)*temp0);
        temptemp = temp(:,1);
        temp = .5*(repmat(temptemp,1,size(val,2)) + temptemp(ind)) - temp;
        D_kernels(:,:,t) = temp+eps;
        t = t+1;
    end
    end
end
end


function   distf = mex_L2_distance_1(F, ind)
[m,n] = size(ind);    
I = repmat([1:m]',1,n);
temp = sum((F(I(:),:) - F(ind(:),:)).^2,2);
distf = zeros(m,n);
distf(:) = temp;
end
