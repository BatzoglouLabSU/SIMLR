"""
    Functions for large scale SIMLR and accuracy checks

    ---------------------------------------------------------------------

    This module contains the following functions:

    save_sparse_csr
    save a sparse csr format input of single-cell RNA-seq data
    load_sparse_csr
    load a sparse csr format input of single-cell RNA-seq data
    nearest_neighbor_search
    Approximate Nearset Neighbor search for every cell
    NE_dn
    Row-normalization of a matrix
    mex_L2_distance
    A fast way to calculate L2 distance
    Cal_distance_memory
    Calculate Kernels in a memory-saving mode
    mex_multipleK
    A fast way to calculate kernels
    Hbeta
    A simple LP method to solve linear weight
    euclidean_proj_simplex
    A fast way to calculate simplex projection
    fast_pca
    A fast randomized pca with sparse input
    fast_minibatch_kmeans
    A fast mini-batch version of k-means
    SIMLR_Large
    A large-scale implementation of our SIMLR
    ---------------------------------------------------------------------

    Copyright 2016 Bo Wang, Stanford University.
    All rights reserved.
    """

from __future__ import print_function
import numpy as np
import sys
import os
from annoy import AnnoyIndex
import scipy.io as sio
from scipy.sparse import csr_matrix, csc_matrix, linalg
from fbpca import svd, pca
import time
from sklearn.decomposition import TruncatedSVD
from sklearn.cluster import MiniBatchKMeans, KMeans
from sklearn.metrics.cluster import normalized_mutual_info_score as nmi
from sklearn.metrics.cluster import adjusted_rand_score as ari


def save_sparse_csr(filename,array, label=[]):
    np.savez(filename,data = array.data ,indices=array.indices,
          indptr =array.indptr, shape=array.shape, label = label )

def load_sparse_csr(filename):
    loader = np.load(filename)
    if 'label' in loader.keys():
        label = loader['label']
    else:
        label = []
    return csr_matrix((  np.log10(1.0+loader['data']), loader['indices'], loader['indptr']),
                      shape = loader['shape']), label

def nearest_neighbor_search(GE_csc, K = 20):
    n,d = GE_csc.shape
    t = AnnoyIndex(d)
    for i in xrange(n):
        t.add_item(i,GE_csc[i,:])
    t.build(2*d)
    t.save('test.ann')
    u = AnnoyIndex(d)
    u.load('test.ann')
    os.remove('test.ann')
    val = np.zeros((n,K))
    ind = np.zeros((n,K))
    for i in xrange(n):
        tmp, tmp1 = u.get_nns_by_item(i,K, include_distances=True)
        ind[i,:] = tmp
        val[i,:] = tmp1
    return ind.astype('int'), val


def NE_dn(A, type='ave'):
    m , n = A.shape
    diags = np.ones(m)
    diags[:] = abs(A).sum(axis=1).flatten()
    if type == 'ave':
        D = 1/(diags+np.finfo(float).eps)
        return A*D[:,np.newaxis]
    elif type == 'gph':
        D = 1/np.sqrt(diags+np.finfo(float).eps)
        return (A*D[:,np.newaxis])*D


def mex_L2_distance(F, ind):
    m,n = ind.shape
    I = np.tile(np.arange(m), n)
    temp = np.take(F, I, axis = 0) - np.take(F, ind.ravel(order = 'F'),axis=0)
    temp = (temp*temp).sum(axis=1)
    return temp.reshape((m,n),order = 'F')



def Cal_distance_memory(S, alpha):
    NT = len(alpha)
    DD = alpha.copy()
    for i in xrange(NT):
        temp = np.load('Kernel_' + str(i)+'.npy')
        if i == 0:
            distX = alpha[0]*temp
        else:
            distX += alpha[i]*temp
        DD[i] = ((temp*S).sum(axis = 0)/(S.shape[0]+0.0)).mean(axis = 0)
    alphaK0 = umkl_bo(DD, 1.0/len(DD));
    alphaK0 = alphaK0/np.sum(alphaK0)
    return distX, alphaK0


def mex_multipleK(val, ind, KK=20, ismemory = 0):
    #val *= val
    m,n=val.shape
    sigma = np.arange(1,2.1,.25)
    allK = (np.arange(np.ceil(KK/2.0), min(n,np.ceil(KK*1.5))+1, np.ceil(KK/10.0))).astype('int')
    if ismemory:
        D_kernels = []
        alphaK = np.ones(len(allK)*len(sigma))/(0.0 + len(allK)*len(sigma))
    else:
        D_kernels = np.zeros((m,n,len(allK)*len(sigma)))
        alphaK = np.ones(D_kernels.shape[2])/(0.0 + D_kernels.shape[2])
    t = 0;
    for k in allK:
        temp = val[:,np.arange(k)].sum(axis=1)/(k+0.0)
        temp0 = .5*(temp[:,np.newaxis] + np.take(temp,ind))
        temp = val/temp0
        temp*=temp
        for s in sigma:
            temp1 =  np.exp(-temp/2.0/s/s)/np.sqrt(2*np.pi)/s/temp0
            #temp1[:] = NE_dn(temp1)
            temptemp = temp1[:, 0]
            temp1[:] = .5*(temptemp[:,np.newaxis] + temptemp[ind]) - temp1
            if ismemory:
                np.save('Kernel_' + str(t), temp1 - temp1.min())
            else:
                D_kernels[:,:,t] = temp1 - temp1.min()
            t = t+1

    return D_kernels, alphaK

def Hbeta(D,beta):
    D = (D-D.min())/(D.max() - D.min() + np.finfo(float).eps)
    P = np.exp(-D*beta)
    sumP = P.sum()
    H = np.log(sumP) + beta*sum(D*P)/sumP
    P = P / sumP
    return H, P



def umkl_bo(D, beta):
    tol = 1e-4
    u = 20
    logU = np.log(u)
    H, P = Hbeta(D,beta)
    betamin = -np.inf
    betamax = np.inf
    Hdiff = H - logU
    tries = 0
    while(abs(Hdiff)>tol)&(tries < 30):
        if Hdiff>0:
            betamin = beta
            if np.isinf(betamax):
                beta *= 2.0
            else:
                beta = .5*(beta + betamax)
        else:
            betamax = beta
            if np.isinf(betamin):
                beta /= 2.0
            else:
                beta = .5*(beta + betamin)
        H, P = Hbeta(D,beta)
        Hdiff = H - logU
        tries +=1
    return P

def euclidean_proj_simplex(v, s=1):
    """ Compute the Euclidean projection on a positive simplex
    Solves the optimisation problem (using the algorithm from [1]):
        min_w 0.5 * || w - v ||_2^2 , s.t. \sum_i w_i = s, w_i >= 0
    Parameters
    ----------
    v: (n,) numpy array,
       n-dimensional vector to project
    s: int, optional, default: 1,
       radius of the simplex
    Returns
    -------
    w: (n,) numpy array,
       Euclidean projection of v on the simplex
    """
    assert s > 0, "Radius s must be strictly positive (%d <= 0)" % s
    n,d = v.shape  # will raise ValueError if v is not 1-D
    # get the array of cumulative sums of a sorted (decreasing) copy of v

    v -= (v.mean(axis = 1)[:,np.newaxis]-1.0/d)
    u = -np.sort(-v)
    cssv = np.cumsum(u,axis = 1)
    # get the number of > 0 components of the optimal solution
    temp = u * np.arange(1,d+1) - cssv +s
    temp[temp<0] = 'nan'
    rho = np.nanargmin(temp,axis = 1)
    #rho = np.nonzero(u * np.arange(1, n+1) > (cssv - s))[0][-1]
    # compute the Lagrange multiplier associated to the simplex constraint
    theta = (cssv[np.arange(n), rho] - s) / (rho + 1.0)
    # compute the projection by thresholding v using theta
    w = (v - theta[:,np.newaxis]).clip(min=0)
    return w

def fast_eigens(val, ind,k):
    n,d = val.shape
    rows = np.tile(np.arange(n), d)
    cols = ind.ravel(order='F')
    A = csr_matrix((val.ravel(order='F'),(rows,cols)),shape = (n, n)) + csr_matrix((val.ravel(order='F'),(cols,rows)),shape = (n, n))
    (d,V) = linalg.eigsh(A,k,which='LM')
    d = -np.sort(-np.real(d))
    return np.real(V),d

def fast_pca(in_X, no_dim):
    (U, s, Va) = pca(in_X, no_dim, True, 5)
    del Va
    U[:] = U*np.sqrt(np.abs(s))
    D = 1/(np.sqrt(np.sum(U*U,axis = 1)+np.finfo(float).eps)+np.finfo(float).eps)
    return U*D[:,np.newaxis]

def fast_minibatch_kmeans(X,C):
    cls = MiniBatchKMeans(n_clusters=C, n_init = 100, max_iter = 100)
    return cls.fit_predict(X)

def SIMLR_Large(X, c, K = 20,  is_memory = 0, NITER = 5, beta = 0.8):
    n,d = X.shape
    if d > 500:
        print('SIMLR highly recommends you to perform PCA first on the data\n');
        print('Please use the in-line function fast_pca on your input\n');
    ind, val = nearest_neighbor_search(X, min(K*2, np.ceil(n/10.0)))
    del X
    D_Kernels, alphaK = mex_multipleK(val, ind, K, is_memory)
    del val
    if is_memory:
        distX,alphaK = Cal_distance_memory(np.ones((ind.shape[0], ind.shape[1])), alphaK)
    else:
        distX = D_Kernels.dot(alphaK)
    rr = (.5*(K*distX[:,K+2] - distX[:,np.arange(1,K+1)].sum(axis = 1))).mean()
    lambdar = rr
    S0 = distX.max() - distX
    S0[:] = NE_dn(S0)
    F, evalues = fast_eigens(S0.copy(), ind.copy(), c)
    #F = NE_dn(F)
    F0 = F.copy()
    for iter in range(NITER):
        FF = mex_L2_distance(F, ind)
        FF[:] = (distX + lambdar*FF)/2.0/rr
        #FF[:] = euclidean_proj_simplex(-FF.transpose()).transpose()
        FF[:] = euclidean_proj_simplex(-FF)
        #FF[:] = NE_dn(FF)
        S0[:] = (1-beta)*S0 + beta*FF

        F[:], evalues = fast_eigens(S0, ind,c)
        #F[:] = NE_dn(F)
        F[:] = (1-beta)*F0 + beta*F
        F0 = F.copy()
        lambdar = lambdar * 1.5
        rr = rr / 1.05
        if is_memory:
            distX, alphaK0 = Cal_distance_memory(S0, alphaK)
            alphaK = (1-beta)*alphaK + beta*alphaK0
            alphaK = alphaK/np.sum(alphaK)

        else:
            DD = ((D_Kernels*S0[:,:,np.newaxis]).sum(axis = 0)/(D_Kernels.shape[0]+0.0)).mean(axis = 0)
            alphaK0 = umkl_bo(DD, 1.0/len(DD));
            alphaK0 = alphaK0/np.sum(alphaK0)
            alphaK = (1-beta)*alphaK + beta*alphaK0
            alphaK = alphaK/np.sum(alphaK)
            distX = D_Kernels.dot(alphaK)

    if is_memory:
        for i in xrange(len(alphaK)):
            os.remove('Kernel_' + str(i) + '.npy')
    rows = np.tile(np.arange(n), S0.shape[1])
    cols = ind.ravel(order='F')
    val = S0
    S0 = csr_matrix((S0.ravel(order='F'),(rows,cols)),shape = (n, n)) + csr_matrix((S0.ravel(order='F'),(cols,rows)),    shape = (n, n))
    return S0, F, val, ind



if __name__ == "__main__":
    filename = 'demo/Zeisel.npz'
    GE, label = load_sparse_csr(filename)
    c = label.max()
    thres = 0.1
    stdsparsity = GE.getnnz(axis = 1)/(0.0 + GE.shape[1])
    X = GE[np.where(stdsparsity > thres)[0],:].T
    del GE
    if X.shape[1]>500:
        X = fast_pca(X,500)
    else:
        X = X.todense()
    print('Preprocess Done! Start to Run SIMLR!\n')
    start_main = time.time()
    S, F, val, ind = SIMLR_Large(X, c, 30)
    print('Successfully Run SIMLR! SIMLR took %f seconds in total\n' % (time.time() - start_main))
    y_pred = fast_minibatch_kmeans(F,c)
    print('NMI value is %f \n' % nmi(y_pred.flatten(),label.flatten()))
    print('ARI value is %f \n' % ari(y_pred.flatten(),label.flatten()))

