function D = L2_distance_K_bo(X,sigma)
% if nargin==1
%     D = L2_distance_1(X,X);
%     maxt = max(D);
%     D(D==0) = 10000;
%     mint = min(D);
%     sigma = ((maxt-mint)/2/log(maxt./mint));
%     sigma = sigma'*sigma;
% end
sigma = 2;
%%X is of size Nxd
X = X';
[N,d] = size(X);
%K = getKern(X',X',sigma);
%K1 = X*X';
t = 0;
while t<10;
    if d>100000
        NN = 2;
        ind = union(1:round(d/NN):d,d);
        temp = randperm(d);
        K = zeros(N);
        for i = 1:(length(ind)-1)
            K = K + affinityMatrix(dist2(X(:,temp(ind(i):ind(i+1)))),15,sigma);
        end
    else
        K = affinityMatrix(dist2(X),15,sigma);
    end
    t = t+1;
end
K = K/t;
%K = (0.4*kernel_norm(K1)+0.6*kernel_norm(K2));
%K = kernel_norm(K);
D = repmat(diag(K),1,N) + repmat(diag(K)',N,1) - 2*K;

function G = kernel_norm(K)
k = 1./sqrt(diag(max(K))+0.0001);
G = K.*(k*k');