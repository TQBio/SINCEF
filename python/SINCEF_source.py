# -*- coding: utf-8 -*-
"""
Created on Mon May 10 20:24:45 2021

@author: tianqi
"""

from scipy.spatial.distance import pdist
from scipy.spatial.distance import squareform
from sklearn import metrics
from sklearn.manifold import SpectralEmbedding
from sklearn.cluster import AgglomerativeClustering
import pandas as pd

###Input: S—— the distance matrix
###output: A—— the affinity matrix
def myKNN(S, k, sigma=1.0):
    N = len(S)
    A = np.zeros((N,N))
    for i in range(N):
        dist_with_index = zip(S[i], range(N))
        dist_with_index = sorted(dist_with_index, key=lambda x:x[0])
        neighbours_id = [dist_with_index[m][1] for m in range(k+1)] 
        for j in neighbours_id: # xj is xi's neighbour
            A[i][j] = np.exp(-S[i][j]/2/sigma/sigma)
            A[j][i] = A[i][j] # mutually
    return A

###Input: aff1,aff2,aff3—— the affinity matrices, corresponding to Cosine, Hamming and Pearson distance; 
###       dim—— the number of embedding dimensions;
###output: dism—— the reconstructed novel distance matrix
def SPC_Embedding(aff1,aff2,aff3,dim):
    fs1 = SpectralEmbedding(n_components=dim, random_state=123, affinity='precomputed').fit_transform(aff1)
    fs2 = SpectralEmbedding(n_components=dim, random_state=123, affinity='precomputed').fit_transform(aff2)
    fs3 = SpectralEmbedding(n_components=dim, random_state=123, affinity='precomputed').fit_transform(aff3)
    fs_dis = pd.concat([pd.DataFrame(fs1),pd.DataFrame(fs2),pd.DataFrame(fs3)],axis=1)
    dism = squareform(pdist(fs_dis,metric='euclidean'))
    return(dism)

###Input: dism—— the distance matrix; C—— the number of clusters
###output: label—— the predicted cell label of hierarchical clustering
def hc_pre(dism,C):
    hc = AgglomerativeClustering(n_clusters=C, affinity='precomputed',linkage='average')
    out = hc.fit_predict(dism)
    return out

###Input: c1—— the predicted cell label; c2—— the true cell label
###output: score—— the ARI score for the predicted cell label
def get_ARI(c1,c2):
    score = metrics.adjusted_rand_score(c1,c2)
    return(score)

###Input: c1—— the predicted cell label; c2—— the true cell label
###output: score—— the NMI score for the predicted cell label
def get_NMI(c1,c2):
    score = metrics.normalized_mutual_info_score(c1,c2)
    return(score)

###Input: L—— the predicted cell label; dism—— the distance matrix
###output: score—— the silhouette coefficient of the predicted cell label
def get_silhouette_score(L,dism):
    score = metrics.silhouette_score(dism,L,metric='precomputed')
    return(score)