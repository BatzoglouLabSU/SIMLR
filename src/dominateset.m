function  [PNN_matrix]  = dominateset(aff_matrix,NR_OF_KNN);

%[res,loc] = maxk(aff_matrix, NR_OF_KNN, 2 );
%inds = repmat((1:size(loc,1))',[1 size(loc,2)]);

[A,B] = sort(aff_matrix,2,'descend');
res = A(:,1:NR_OF_KNN);
inds = repmat([1:length(aff_matrix)]',1,NR_OF_KNN);
loc = B(:,1:NR_OF_KNN);
PNN_matrix1 = zeros(size(aff_matrix));
PNN_matrix1(sub2ind(size(aff_matrix),inds(:),loc(:))) = res(:);
PNN_matrix = full(PNN_matrix1+PNN_matrix1')/2;
%PNN_matrix = max(PNN_matrix1,PNN_matrix1');
% 
% N = length(aff_matrix);
% PNN_matrix = bmatch(aff_matrix,NR_OF_KNN*ones(N,1));
% PNN_matrix = aff_matrix.*PNN_matrix;
% PNN_matrix = full(PNN_matrix+PNN_matrix')/2;
% a = sum(PNN_matrix);
% if sum(a<eps)
%     PNN_matrix
% end


end

