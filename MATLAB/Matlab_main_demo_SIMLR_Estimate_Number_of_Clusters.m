clear
clc
close all

addpath('data')
addpath('src')
dataset = {'4_Usoskin'}

load(['Test_' dataset{1}]);
rng(95584,'twister'); %%% for reproducibility
[K1, K2] = Estimate_Number_of_Clusters_SIMLR(in_X,2:10);
