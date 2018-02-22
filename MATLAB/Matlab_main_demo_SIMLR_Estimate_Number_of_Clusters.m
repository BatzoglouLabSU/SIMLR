clear
clc
close all

addpath('data')
addpath('src')
dataset = {'3_Pollen'}

load(['Test_' dataset{1}]);
rng('default'); %%% for reproducibility
[K1, K2] = Estimate_Number_of_Clusters_SIMLR(in_X,2:5);
