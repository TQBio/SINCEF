# SINCEF
Clustering single-cell methylation data using a novel distance based on spectral embedding aggregation. SINCEF contains both R and python scripts to achieve the clustering process.

# Tutorial

## Step 1 

Prepare your single-cell methylation dataset in the specified format. The recommended format is an R list, where each element is a dataframe containing methylation information for a single cell, organized as follows:
   
    /Chrom/  /location/   /methylation_state/
   
    chr1      131009             1
    
    chr1      131031             0
    
    ...        ...              ...
    
    chrY     56770241            1
    
   This repositiry provides the simulated single-cell methylation data ("data_bd.RData") which is sampled from the GSE87197 and more information on the data format can be          obtained from the simulated dataset . 
    
   Download link: https://github.com/TQBio/SINCEF/blob/master/R/data_bd.RData

## Step 2

   Calculate the cell-to-cell distance matrices by parallel computing in R environment.
    
   Example:
    
    load('data_bd.RData') ###Load simulated data
    
    ls()
       
    dism_pearson <- Output_DISM(k_cpu=8,data_BLOOD,method='Pearson')  #k_cpu means the number of processors to be used
        
    dism_cosine  <- Output_DISM(k_cpu=8,data_BLOOD,method='Cosine')
    
    dism_hamming <- Output_DISM(k_cpu=8,data_BLOOD,method='Hamming')
    
    ###Save the resulting distance matrices
    
    write.csv(dism_pearson, file='dism_pearson.csv')
    
    write.csv(dism_cosine, file='dism_cosine.csv')
    
    write.csv(dism_hamming, file='dism_hamming.csv')
    
## Step 3

   Using SINCEF to obtain the reconstructed novel distance matrix and its corresponding hierarchical clustering object.
   
    ###Execute the Python scripts or load all the functions in python compiler environment.
    
    clus = find_clus(dism, k_max=11) ### find the optimal number of clusters c
    
    aff1 = myKNN(dism_cosine, 5) ###calculate the affinity matrix with neighbors = 5
    aff2 = myKNN(dism_hamming, 5)
    aff3 = myKNN(dism_pearson, 5)
    
    rdm = SPC_Embedding(aff1, aff2, aff3,2) ### default embedded dimension = 2;
    
    hc_pre = hc_pre(rdm,c)     ### get the predicted cell labels
    
           
## Step 4

   Evaluate clustering performance using ARI and NMI.
   
    ari_score = get_ARI(hc_pre, true_label)
    
    nmi_score = get_ARI(hc_pre, true_label)
    
    sil_score = get_silhouette_score(hc_pre, rdm)

# Note

  1. When calculating the dissimilarity matrix, please set it reasonably according to the number of available cores of your PC or server, otherwise it may cause memory overflow.
  
  2. Any problems or bugs encountered during the use of SINCEF, please contact tqglowing@std.uestc.edu.cn.
