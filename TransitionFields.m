function [ Wnew ] = TransitionFields( W )
% W = W-diag(diag(W));
% W = (W+W')/2;
zeroindex = find(sum(W,2)==0);
%W = dn(W,'ave');
%W = (W+W')/2+eye(length(W));

% maxW = max(max(W));
% W = (W+W')/2+maxW*eye(length(W));

W = dn(W,'ave');
w = sqrt(sum(W)+eps);
W = W./repmat(w,length(W),1);
W = W*W';
Wnew = W;
Wnew(zeroindex,:) = 0;
Wnew(:,zeroindex) = 0;
end

