function  PNN_matrix = dominateset(aff_matrix,NR_OF_KNN);
% P=aff_matrix;K = NR_OF_KNN;
% [YW,IW1] = sort(P,2,'descend');
% clear YW;
% [n,m]=size(P);
% DP=zeros(size(P));
% temp=repmat((1:n)',1,K);
% I1=(IW1(:,1:K)-1)*m+temp;
% DP(I1(:))=P(I1(:));
% % 
% % [YW,IW1] = sort(P',2,'descend');
% % clear YW;
% % [n,m]=size(P);
% % temp=repmat((1:n)',1,K);
% % I1=(IW1(:,1:K)-1)*m+temp;
% % DP(I1(:))=P(I1(:));
% % %DP=DP./repmat(sum(DP,2),1,n);
%  PNN_matrix = DP;



%[res,loc] = maxk(aff_matrix, NR_OF_KNN, 2 );
%inds = repmat((1:size(loc,1))',[1 size(loc,2)]);

%PNN_matrix1 = sparse(inds(:),loc(:),res(:),size(aff_matrix,1), ...
%    size(aff_matrix,2),numel(res));

[A,B] = sort(aff_matrix,2,'descend');
res = A(:,1:NR_OF_KNN);
inds = repmat([1:length(aff_matrix)]',1,NR_OF_KNN);
loc = B(:,1:NR_OF_KNN);
PNN_matrix1 = zeros(size(aff_matrix));
PNN_matrix1(sub2ind(size(aff_matrix),inds(:),loc(:))) = res(:);
PNN_matrix = full(PNN_matrix1+PNN_matrix1')/2;


PNN_matrix = full(PNN_matrix1+PNN_matrix1')/2;

end

