function [K1, K2] = Estimate_Number_of_Clusters_SIMLR(X, NUMC)

D_Kernels = multipleK(X);
distX = mean(D_Kernels,3);
W =  max(max(distX)) - distX;
W = Network_Diffusion(W,max(ceil(size(X,1)/20),10));

[Quality] = Estimate_Number_of_Clusters_given_graph(W, NUMC);
[Quality_plus] = Estimate_Number_of_Clusters_given_graph(W, NUMC+1);
[Quality_minus] = Estimate_Number_of_Clusters_given_graph(W, NUMC-1);

K1 = 2*(1 + Quality) - (2 + Quality_plus + Quality_minus);
K2 = K1.*(NUMC+1)./(NUMC);
subplot(1,2,1)
plot(NUMC,K1,'b-s','LineWidth',4);
title('Relative Quality')
subplot(1,2,2)
plot(NUMC,K2,'r-o','LineWidth',4);
title('Adjusted Quality')

end

function [quality] = Estimate_Number_of_Clusters_given_graph(W, NUMC)
%%%This function estimates the number of clusters given the two huristics
%%%given in the supplementary materials of our nature method paper
%W is the similarity graph
%NUMC is a vector which contains the possible choices of number of
%clusters.

if nargin < 2
    NUMC = 2:5;
end

if min(NUMC)==1
    warning('Note that we always assume there are more than one cluster.');
    NUMC(NUMC<1) = [];
end
W = (W + W')/2;

if ~isempty(NUMC)
    degs = sum(W, 2);
    D    = sparse(1:size(W, 1), 1:size(W, 2), degs);
    % compute unnormalized Laplacian
    L = D - W;
    degs(degs == 0) = eps;
    % calculate D^(-1/2)
    D = spdiags(1./(degs.^0.5), 0, size(D, 1), size(D, 2));
    % calculate normalized Laplacian
    L = D * L * D;
    % compute the eigenvectors corresponding to the k smallest
    % eigenvalues
    [U, eigenvalue] = eig(L);
    eigenvalue  = diag(eigenvalue);
    [a,b] = sort((eigenvalue),'ascend');
    eigenvalue = (eigenvalue(b));
    U = U(:,b);
    eigengap = abs(diff(eigenvalue));
    for ck = NUMC
        Cindex = find(NUMC==ck);
        if ck == 1
            quality(Cindex) = sum(sum(diag(1./(U(:,1)+eps))*U(:,1)));
        else
            UU = U(:,1:ck);
            UU = UU./repmat(sqrt(sum(UU.^2,2))+eps,1,size(UU,2));
            [EigenvectorsDiscrete,EigenVectors ]=discretisation(UU);
            EigenVectors = EigenVectors.^2;
            [temp1,temp] = sort(EigenVectors,2, 'descend');
            quality(Cindex) = (1-eigenvalue(ck+1))/(1-eigenvalue(ck))*sum(sum(diag(1./(temp1(:,1)+eps))*temp1(:,1:max(2,ck-1))));
        end
    end
    
end


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
    k = 1./sqrt(diag(K)+eps);
    G = K.*(k*k');
    
    D_Kernels(:,:,i) = (repmat(diag(G),1,length(G)) +repmat(diag(G)',length(G),1) - 2*G)/2;
    D_Kernels(:,:,i) = D_Kernels(:,:,i) - diag(diag(D_Kernels(:,:,i)));
    
end

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
alpha = 0.8;
beta = 2;
d = (1-alpha)*d./(1-alpha*d.^beta);
D = diag(real(d));
W = U*D*U';
W = (W.*(1-eye(length(W))))./repmat(1-diag(W),1,length(W));
D=sparse(1:length(DD),1:length(DD),DD);

W = (W+W')/2;

end


function [EigenvectorsDiscrete,EigenVectors]=discretisation(EigenVectors)
%

[n,k]=size(EigenVectors);

vm = sqrt(sum(EigenVectors.*EigenVectors,2));
EigenVectors = EigenVectors./repmat(vm+eps,1,k);

R=zeros(k);
R(:,1)=EigenVectors(round(n/2),:)';

c=zeros(n,1);
for j=2:k
    c=c+abs(EigenVectors*R(:,j-1));
    [minimum,i]=min(c);
    R(:,j)=EigenVectors(i,:)';
end

lastObjectiveValue=0;
exitLoop=0;
nbIterationsDiscretisation = 0;
nbIterationsDiscretisationMax = 20;
while exitLoop== 0
    nbIterationsDiscretisation = nbIterationsDiscretisation + 1 ;
    EigenvectorsDiscrete = discretisationEigenVectorData(EigenVectors*R);
    [U,S,V] = svd(EigenvectorsDiscrete'*EigenVectors+eps,0);
    NcutValue=2*(n-trace(S));
    
    if abs(NcutValue-lastObjectiveValue) < eps | nbIterationsDiscretisation > nbIterationsDiscretisationMax
        exitLoop=1;
    else
        lastObjectiveValue = NcutValue;
        R=V*U';
    end
end
end

function Y = discretisationEigenVectorData(EigenVector)
% Y = discretisationEigenVectorData(EigenVector)
%
% discretizes previously rotated eigenvectors in discretisation
% Timothee Cour, Stella Yu, Jianbo Shi, 2004

[n,k]=size(EigenVector);


[Maximum,J]=max(EigenVector');

Y=sparse(1:n,J',1,n,k);

end
