% This script is used to compile a few C-mex files necessary to run SIMLR large scale 

% first of all, mex projsplx_c.c and Kbeta.cpp 
cd src/ 
mex -largeArrayDims projsplx_c.c 
mex -largeArrayDims Kbeta.cpp 

% second mex tSNE 
mex -largeArrayDims compute_wtsne_obj_grad_repulsive_barneshut.c barnes_hut.c 

% next we need to compile the eigen-solver with the external library Spectra 
mex -largeArrayDims mex_top_eig.cpp -I../External/spectra-master/include/ 

% finally we need to compile Annoy library for ANN search 
mex CXXFLAGS="\$CXXFLAGS -std=c++11 -lm -lgsl -lgslcblas -O3 -ffast-math" mex_KNN_Annoy.cpp -largeArrayDims -I../External/ 

% Note that if you have trouble compiling the above comand by missing a few 
% library, please use the following one instead by uncommenting (this will potential slow the ANN search) 
% mex mex_KNN_Annoy.cpp -largeArrayDims -I../External/ 

% you have succesfully compile all the files. Exit to the main folder 
cd ../
addpath('src')
