function [y, S, F, ydata,alphaK,timeOurs,converge,LF] = CIMLR(alldata, c, k)

if nargin==2
    k=10;
    ifimpute = 0;
    normalize = 0;
end



t0 = tic;
order = 2;
no_dim=c;


NITER = 30;
num = size(alldata{1},1);
r = -1;
beta = 0.8;
for i = 1:length(alldata)
    
    if i == 1
        D_Kernels = multipleK(alldata{i});
    else
        D_Kernels = cat(3, D_Kernels,multipleK(alldata{i}));
    end
end
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
S0 = max(max(distX))-distX;
S0 = Network_Diffusion(S0,k);
S0 = NE_dn(S0,'ave');
S= (S0 + S0')/2;
D0 = diag(sum(S,order));
L0= D0-S;
[F, temp, evs]=eig1(L0, c, 0);
F = NE_dn(F,'ave');
for iter = 1:NITER
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
    S = (1-beta)*A+beta*S;
    S = Network_Diffusion(S,k);
    S= (S + S')/2;
    D = diag(sum(S,order));
    L = D - S;
    F_old = F;
    [F, temp, ev]=eig1(L, c, 0);
    F = NE_dn(F,'ave');
    F = (1-beta)*F_old+beta*F;
    evs(:,iter+1) = ev;
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
    converge(iter) = fn2-fn1;
    if iter<10
        if (ev(end) > 0.000001)
            lambda = 1.5*lambda;
            r = r/1.01;
        end
    else
        if (converge(iter)>1.01*converge(iter-1))
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
end

LF = F;
D = diag(sum(S,order));
L = D - S;
[U,D] = eig(L);
if length(no_dim)==1
    F = tsne_p_bo((S),[], U(:,1:no_dim));
else
    clear F;
    for i = 1:length(no_dim)
        F{i} = tsne_p_bo((S),[], U(:,1:no_dim(i)));
    end
end
timeOurs = toc(t0);

[~,center] = litekmeans(LF, c,'replicates',200);
[~,center] = min(dist2(center,LF),[],2);
y = litekmeans(F,c,'Start',center);
ydata = tsne_p_bo(S);

end

function thisP = umkl_bo(D,beta)
if nargin<2
    beta = 1/length(D);
end
tol = 1e-4;
u = 150;logU = log(u);
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


N = size(x,1);
KK = 0;
sigma = [2:-0.25:1];
Diff = (dist2(x));
[T,INDEX]=sort(Diff,2);
[m,n]=size(Diff);
allk = 10:2:30;
t=1;
for l = 1:length(allk)
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

for i = 1:size(Kernels,3)
    K = Kernels(:,:,i);
    k = 1./sqrt(diag(K)+1);
    G = K.*(k*k');
    %G = K;
    D_Kernels(:,:,i) = (repmat(diag(G),1,length(G)) +repmat(diag(G)',length(G),1) - 2*G)/2;
    D_Kernels(:,:,i) = D_Kernels(:,:,i) - diag(diag(D_Kernels(:,:,i)));
  
end

end



%
function W = Network_Diffusion(A, K)
%K = min(2*K, round(length(A)/10));
A = A-diag(diag(A));
P = (dominateset(double(abs(A)),min(K,length(A)-1))).*sign(A);
DD = sum(abs(P'));
P = P + (eye(length(P))+diag(sum(abs(P'))));
P = (TransitionFields(P));
[U,D] = eig(P);
d = real((diag(D))+eps);
alpha = 0.8;
beta = 2;
d = (1-alpha)*d./(1-alpha*d.^beta);


D = diag(real(d));
W = U*D*U';

W = (W.*(1-eye(length(W))))./repmat(1-diag(W),1,length(W));
D=sparse(1:length(DD),1:length(DD),DD);
W=D*(W);
W = (W+W')/2;

end

function ydata = tsne_p_bo(P, labels, no_dims)
%TSNE_P Performs symmetric t-SNE on affinity matrix P
    if ~exist('labels', 'var')
        labels = [];
    end
    if ~exist('no_dims', 'var') || isempty(no_dims)
        no_dims = 2;
    end
    
    % First check whether we already have an initial solution
    if numel(no_dims) > 1
        initial_solution = true;
        ydata = no_dims;
        no_dims = size(ydata, 2);
    else
        initial_solution = false;
    end
    
    % Initialize some variables
    n = size(P, 1);                                     % number of instances
    momentum = .08;                                     % initial momentum
    final_momentum = .1;                               % value to which momentum is changed
    mom_switch_iter = 250;                              % iteration at which momentum is changed
    stop_lying_iter = 100;                              % iteration at which lying about P-values is stopped
    max_iter = 1000;                                    % maximum number of iterations
    epsilon = 500;                                      % initial learning rate
    min_gain = .01;                                     % minimum gain for delta-bar-delta
    
    % Make sure P-vals are set properly
    P(1:n + 1:end) = 0;                                 % set diagonal to zero
    P = 0.5 * (P + P');                                 % symmetrize P-values
    P = max(P ./ sum(P(:)), realmin);                   % make sure P-values sum to one
    const = sum(P(:) .* log(P(:)));                     % constant in KL divergence
    if ~initial_solution
        P = P * 4;                                      % lie about the P-vals to find better local minima
    end
    
    % Initialize the solution
    if ~initial_solution
        ydata = .0001 * randn(n, no_dims);
    end
    y_incs  = zeros(size(ydata));
    gains = ones(size(ydata));
    
    % Run the iterations
    for iter=1:max_iter
        
        % Compute joint probability that point i and j are neighbors
        sum_ydata = sum(ydata .^ 2, 2);
        num = 1 ./ (1 + bsxfun(@plus, sum_ydata, bsxfun(@plus, sum_ydata', -2 * (ydata * ydata')))); % Student-t distribution
        num(1:n+1:end) = 0;                                                 % set diagonal to zero
        Q = max(num ./ sum(num(:)), realmin);                               % normalize to get probabilities
        
        % Compute the gradients (faster implementation)
        L = (P - Q) .* num;
        y_grads = 4 * (diag(sum(L, 1)) - L) * ydata;
            
        % Update the solution
        gains = (gains + .2) .* (sign(y_grads) ~= sign(y_incs)) ...         % note that the y_grads are actually -y_grads
              + (gains * .8) .* (sign(y_grads) == sign(y_incs));
        gains(gains < min_gain) = min_gain;
        y_incs = momentum * y_incs - epsilon * (gains .* y_grads);
        ydata = ydata + y_incs;
        ydata = bsxfun(@minus, ydata, mean(ydata, 1));
        ydata(ydata<-100)=-100;ydata(ydata>100)=100;
        
        % Update the momentum if necessary
        if iter == mom_switch_iter
            momentum = final_momentum;
        end
        if iter == stop_lying_iter && ~initial_solution
            P = P ./ 4;
        end
        
        % Print out progress
        if ~rem(iter, 10)
            cost = const - sum(P(:) .* log(Q(:)));
            disp(['Iteration ' num2str(iter) ': error is ' num2str(cost)]);
        end
        
        % Display scatter plot (maximally first three dimensions)
        if ~isempty(labels)
            if no_dims == 1
                scatter(ydata, ydata, 9, labels, 'filled');
            elseif no_dims == 2
                scatter(ydata(:,1), ydata(:,2), 9, labels, 'filled');
            else
                scatter3(ydata(:,1), ydata(:,2), ydata(:,3), 40, labels, 'filled');
            end
            axis equal tight
%             axis off
            drawnow
        end
    end
end
