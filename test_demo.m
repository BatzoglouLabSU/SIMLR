clc
clear 
close all
load mECS_data
[y, A, F, ydata,evs,timeOurs,converge] = GB_OURS_fast((X)', max(label));
%[cluster_result, classification_result] = GB_get_metric(F,label);

%%
%load Kolo_data
%[y, A, F, ydata,evs,timeOurs,converge] = GB_OURS((X)', max(label));
%[cluster_result, classification_result] = GB_get_metric(F,label);

