clear
clc
close all

addpath('data')
addpath('src')
dataset = {'gliomas_multi_omic_data'}

% import the data
load(['Test_6_' dataset{1}]);
tumor_data = gliomas_multi_omic_data;
alldata{1} = tumor_data.point_mutations;
cna_linear = tumor_data.cna_linear;
cna_linear(cna_linear>10) = 10;
cna_linear(cna_linear<-10) = -10;
alldata{2} = cna_linear;
alldata{3} = tumor_data.methylation;
expression = tumor_data.expression;
expression(expression>10) = 10;
expression(expression<-10) = -10;
alldata{4} = tumor_data.expression;
for i = 1:size(alldata{1},2)
    alldata{1}(:,i) = (alldata{1}(:,i) - min(alldata{1}(:,i))) / (max(alldata{1}(:,i)) - min(alldata{1}(:,i)));
    alldata{1}(isnan(alldata{1}(:,i)),i) = 0.5;
    alldata{2}(:,i) = (alldata{2}(:,i) - min(alldata{2}(:,i))) / (max(alldata{2}(:,i)) - min(alldata{2}(:,i)));
    alldata{2}(isnan(alldata{2}(:,i)),i) = 0.5;
    alldata{3}(:,i) = (alldata{3}(:,i) - min(alldata{3}(:,i))) / (max(alldata{3}(:,i)) - min(alldata{3}(:,i)));
    alldata{3}(isnan(alldata{3}(:,i)),i) = 0.5;
    alldata{4}(:,i) = (alldata{4}(:,i) - min(alldata{4}(:,i))) / (max(alldata{4}(:,i)) - min(alldata{4}(:,i)));
    alldata{4}(isnan(alldata{4}(:,i)),i) = 0.5;
end

% estimate the best number of clusters
NUMC = 2:15;
rng(32655,'twister'); %%% for reproducibility
[K1, K2] = Estimate_Number_of_Clusters_CIMLR(alldata,NUMC);
