function [ Wnew ] = TransitionFields( W )
zeroindex = find(sum(W,2)==0);
W = W*length(W);
W = NE_dn(W,'ave');
w = sqrt(sum(abs(W))+eps);
W = W./repmat(w,length(W),1);
W = W*W';
Wnew = W;
Wnew(zeroindex,:) = 0;
Wnew(:,zeroindex) = 0;
end

