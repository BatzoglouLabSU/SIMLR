function Y = SIMLR_embedding_tsne(P, do_init,DD, Y0)
%%%compute 2-D coordinates by approximae t-SNE
%  
%   Y = SIMLR_embedding_tsne(P, do_init)
% 
% Input: 
%   P: N x N, pairwise (sparse) similarities or network (weighted) adjacency matrix
%   do_init: boolean, do over-attraction initialization or not, default=false
% 
% Output:
%   Y: N x 2, 2-D coordinates
%
% All rights reserved.

if ~exist('do_init', 'var') || isempty(do_init)
    do_init = false;
end
if nargin<3
    DD = 2;
end
method = 'wtsne';
P = 0.5 * (P+P');
P = P / sum(sum(P));
theta = 2;
rng('default');
if nargin<4
    Y0 = randn(size(P,1),DD)*1e-4;
end
check_step = 1;
tol = 1e-4;
max_time = inf;
verbose = false;
optimizer = 'fphssne';
recordYs = false;

if do_init
    fprintf('===============================================\n');
    fprintf('Starting Learning tSNE embeddings\n');
    fprintf('===============================================\n');
    attr = 4;
    max_iter = 30;
    Y1 = tSNE_embed(P, Y0, attr, theta, max_iter, check_step, tol, max_time, verbose, recordYs);
else
    Y1 = Y0;
end

attr = 1;
max_iter = 300;
Y = tSNE_embed(P, Y1, attr, theta, max_iter, check_step, tol, max_time, verbose, recordYs);
end
function [Y, ts, Ys] = tSNE_embed(P, Y0, attr, theta, max_iter, check_step, tol, max_time, verbose, recordYs)

Y = Y0;

[I,J] = find(P);
Pnz = nonzeros(P);
weights = sum(P)';

ts(1) = 0;
t = 1;
if recordYs
    Ys = cell(max_iter,1);
    Ys{1} = Y;
else
    Ys = [];
end
tic;
n = size(P,1);

for iter=2:max_iter
    Y_old = Y;
    qnz = 1./(1+sum((Y(I,:)-Y(J,:)).^2,2));
    Pq = attr * sparse(I,J,Pnz.*qnz,n,n);
    [~, repu] = compute_wtsne_obj_grad_repulsive_barneshut(Y, weights, theta, 2);
    
    Y = bsxfun(@rdivide, Pq * Y - repu/4, sum(Pq,2)+eps);
    
    
    if recordYs && mod(iter, check_step)==0
        t = t + 1;
        Ys{t} = Y;
    end
    
    ts(iter) = toc;
    
    if ts(iter)>max_time
        if verbose
            fprintf('max_time exceeded\n');
        end
        break;
    end
    
    diffY = norm(Y-Y_old,'fro') / norm(Y_old,'fro');
    if verbose && mod(iter, check_step)==0
        fprintf('iter=%d/%d, collapsed_time=%.2f, diffY=%.12f, obj=%.12f\n', iter, max_iter, toc, diffY, obj);
    end
    
    if diffY<tol
        if verbose
            fprintf('converged after %d iterations.\n', iter);
        end
        break;
    end
end

ts = ts(1:iter);
if recordYs
    Ys = Ys(1:t);
end
end
