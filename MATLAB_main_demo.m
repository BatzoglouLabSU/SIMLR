clear
clc
close all

addpath(genpath('../'))
dataset = {'1_mECS', '2_Kolod', '3_Pollen', '4_Usoskin'} %% four datasets tested on the paper

for i = 1:4
    load(['Test_' dataset{i}]);
    C = max(true_labs); %%% number of clusters
    rng('default'); %%% for reproducibility
    [y, S, F, ydata] = SIMLR(in_X',C);
    
    %%% report NMI values
    NMI_i = Cal_NMI(y,true_labs);
    fprintf(['The NMI value for dataset ' dataset{i} ' is %f\n'], NMI_i);
    
    %%% visualization
    figure;
    gscatter(ydata(:,1),ydata(:,2),true_labs);
end
