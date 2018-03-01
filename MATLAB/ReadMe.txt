We here provide the Matlab implementations of both SIMLR and CIMLR. 

**SETUP**

Before running our tool, the user needs to mex compile serveral C-mex files. This can be done by running the code SETUP.m. 

Note: the compiling steps implemented in the file SETUP.m are required ONLY to run SIMLR_Large_Scale. All the other scripts do not need them. In order to run the mex command (and, hence, the script SETUP.m), one should have installed a valid compiler. 

**DEMOS**

We provide 3 demos for the usage of SIMLR in both small scale and large scale implementations. In Matlab_main_demo_SIMLR.m we run SIMLR on 4 small-scale datasets while in Matlab_main_demo_SIMLR_Large_Scale.m, we show an example of how to run our tool on 1 large scale dataset. The script Matlab_main_demo_SIMLR_Estimate_Number_of_Clusters.m provides an example of estimation of the best number of clusters by SIMLR. 

Furthermore, we also provide 2 demos for CIMLR. Specifically, the scripts Matlab_main_demo_CIMLR_Estimate_Number_of_Clusters.m and Matlab_main_demo_CIMLR.m show respectively the estimation of the best number of clusters and the subsequent analysis for this number by CIMLR. 
