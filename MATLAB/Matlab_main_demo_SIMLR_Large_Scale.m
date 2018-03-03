clear
clc
close all

addpath('data')
addpath('src')
dataset = {'5_Zeisel'}

% perform the analysis
load(['Test_' dataset{1}]);
C = max(true_labs); %%% number of clusters
rng(75584,'twister'); %%% for reproducibility
in_X = log10(1 + in_X);%%% first of all, take log10 transformation for gene counts
tic
in_X = fast_pca(in_X,500); %%% second, take randomized pca to speed up the ANN search 
toc
tic
[S,F] = SIMLR_Large_Scale(in_X,C,30); %%% S is the learned similarity, F is the latent embedding 
toc
tic
y = litekmeans(F,C,'Replicates',50); %%% run k-means on embeddings to get cell populations
toc

% report NMI values
NMI_i = Cal_NMI(y,true_labs);
fprintf(['The NMI value for dataset ' dataset{1} ' is %f\n'], NMI_i);

% visualization
tic
ydata = SIMLR_embedding_tsne(S,1,2,F(:,1:2));
toc
SIMLR_DisplayVisualization(ydata,true_labs); %%% show visulization
