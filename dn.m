function wn=dn(w,type);
% function Wn=dn(W,type);
% normalizes a symmetric kernel.  
% D=sum(W).  'ave' returns D^-1*W
% 'graph' returns  D^-1/2*W*D^-1/2
% beltrami
D=sum(w,2);

if type == 'ave'
    D=1./D;
    D=sparse(1:length(D),1:length(D),D);
    wn=D*w;
elseif type == 'gph'
    D=1./sqrt(D);
    D=sparse(1:length(D),1:length(D),D);
    wn=D*(w*D);
end