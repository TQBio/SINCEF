#'@title process function
#'
#'@description Obtain common methylation sites for any two cells
#'
#'@param dataframe
#'
#'@return dataframe
#'
#'@examples df_test <- mat_F(df)
#'
#'@export
mat_F <- function(df){
  chr_P <- paste(df[,1], df[,2], sep = "_")
  df <- data.frame(chr_P,df)
  df <- df[,c(1,4)]
  return(df)
  rm(chr_P)
}

#'@title Process function
#'
#'@description Get the corresponding dissimilarity matrix using single-cell methylation profiles
#'
#'@param list
#'
#'@return Pearson similarity matrix
Get_Pearson <- function(j,data){

  cor_cp <- function(i,j,list){
    dm1 <- mat_F(list[[j]])
    dm2 <- mat_F(list[[i]])
    dm <- merge(dm1,dm2,by="chr_P",all=F)
    pv <- cor(dm[,2],dm[,3])
    return(pv)
    rm(dm1)
    rm(dm2)
    rm(dm)
  }###calculate Pearson correlation coefficient

  res_cor <- lapply(1:length(data),cor_cp,j,data)
  res_cor <- as.data.frame(unlist(res_cor))
  return(res_cor)
}

#'@title Process function
#'
#'@description Get the corresponding similarity matrix using single-cell methylation profiles
#'
#'@param list
#'
#'@return Cosine dissimilarity object
Get_Cosine <- function(j,data){

  cor_cos <- function(i,j,list){
    dm1 <- mat_F(list[[j]])
    dm2 <- mat_F(list[[i]])
    dm <- merge(dm1,dm2,by="chr_P",all=F)
    cos <- sum(dm[,2]*dm[,3])/sqrt((sum(dm[,2]^2)*sum(dm[,3]^2)))
    return(cos)
    rm(dm1)
    rm(dm2)
    rm(dm)
  }###calculate cosine correlation coefficient

  res_cor <- lapply(1:length(data),cor_cos,j,data)
  res_cor <- as.data.frame(unlist(res_cor))
  return(res_cor)
}

#'@title Process function
#'
#'@description Get the corresponding dissimilarity matrix using single-cell methylation data
#'
#'@param list
#'
#'@return Dual dissimilarity object
Get_Hamming <- function(j,data){

  cor_ham <- function(i,j,list){
    dm1 <- mat_F(list[[j]])
    dm2 <- mat_F(list[[i]])
    dm <- merge(dm1,dm2,by="chr_P",all=F)
    dm$cot <- ifelse(dm[,2]==dm[,3],1,0)
    dua <- sum(dm$cot)/length(dm$cot)
    return(dua)
    rm(dm1)
    rm(dm2)
    rm(dm)
  }###calculate cosine correlation coefficient

  res_cor <- lapply(1:length(data),cor_ham,j,data)
  res_cor <- as.data.frame(unlist(res_cor))
  return(res_cor)
}

#'@title similarity object
#'
#'@description Get the corresponding dissimilarity matrix using single-cell methylation data
#'
#'@param k_cpu number of CPU cores for parallel computing
#'
#'@param data R list which is composed of dataframes of cell methylation information, including "Chromosome number", "locus position", and "methylation_level(0 or 1)"
#'
#'@param method the distance measure to be used. This must be one of "Cosine" ,"Dual" or "Pearson"
#'
#'@return the corresponding dissimilarity matrix
#'
#'@examples
#'data_cell <- load('data_blood.RData')
#'view(data_cell)
#'dism_pearson <- Output_DISM(k_cpu=8,data_cell,method='Pearson')
#'
#'@export
#'
#'@import parallel
Output_DISM <- function(k_cpu,data,method){
  if(method=='Cosine'){
    cl <- makeCluster(k_cpu)
    clusterExport(cl,"mat_F",envir = environment())
    results <- parLapply(cl,1:length(data),Get_Cosine,data)
    stopCluster(cl)
    res_mat <- do.call(cbind,results)
    res_mat <- matrix(1,nrow = length(data),ncol = length(data))-res_mat
    names(res_mat) <- c(1:length(data))
    return(res_mat)
  }
  if(method=='Hamming'){
    cl <- makeCluster(k_cpu)
    clusterExport(cl,"mat_F",envir = environment())
    results <- parLapply(cl,1:length(data),Get_Dual,data)
    stopCluster(cl)
    res_mat<-do.call(cbind,results)
    res_mat <- matrix(1,nrow = length(data),ncol = length(data))-res_mat
    names(res_mat) <- c(1:length(data))
    return(res_mat)
  }
  if(method=='Pearson'){
    cl <- makeCluster(k_cpu)
    clusterExport(cl,"mat_F",envir = environment())
    results <- parLapply(cl,1:length(data),Get_Pearson,data)
    stopCluster(cl)
    res_mat<-do.call(cbind,results)
    res_mat <- matrix(1,nrow = length(data),ncol = length(data))-res_mat
    names(res_mat) <- c(1:length(data))
    return(res_mat)
  }
  else{
    print('Something wrong with your settings,please check...')
  }

}

#'@title Dissimilarity matrix reconstruction
#'
#'@description The input distance matrices are embedded and fused, and ouput the normalized reconstructed distance matrix and its corresponding hierarchical clustering object.
#'
#'@param d1 one of the three dissimilarity matrices, including Pearson_dism, Cosine_dism and Dual_dism
#'
#'@param d2 one of the three dissimilarity matrices, including Pearson_dism, Cosine_dism and Dual_dism
#'
#'@param d3 one of the three dissimilarity matrices, including Pearson_dism, Cosine_dism and Dual_dism
#'
#'@param is_scale whether to use normalization for reconstructing dissimilarity matrix.
#'
#'@param dim embedded dimensions for initial dissimilarity matrices, recommended between 2 and 5.
#'
#'@return list which is composed of reconstructed dissimilarity matrix,and its hierarchical clustering object.
#'
#'@examples
#'
#'dism_pearson <- Output_DISM(k_cpu=8,data_cell,method='Pearson')
#'dism_cosine  <- Output_DISM(k_cpu=8,data_cell,method='Cosine')
#'dism_dual    <- Output_DISM(k_cpu=8,data_cell,method='Dual')
#'rsm_obj      <- Get_resm(dism_cosine,dism_dual,dism_pearson,is_scale=T)
#'
#'@export
#'
#'@import uwot
Get_resm <- function(d1,d2,d3,dim,is_scale){
  f1 <- umap(as.dist(d1),n_components = dim)
  f2 <- umap(as.dist(d2),n_components = dim)
  f3 <- umap(as.dist(d3),n_components = dim)
  if(is_scale==TRUE){
    f1 = nm(f1)
    f2 = nm(f2)
    f3 = nm(f3)
  }
  else{
    f1 = f1
    f2 = f2
    f3 = f3
  }
  dm1<- as.matrix(dist(f1))
  dm2<- as.matrix(dist(f2))
  dm3<- as.matrix(dist(f3))
  res <- 1/3 * (dm1+dm2+dm3)
  hc_obj <- hclust(as.dist(res),method = 'ward.D2')
  out <- list(res,hc_obj)
  names(out) <- c("RCDM","HC_obj")
  return(out)
}

#'@title Hierarchical clustering visualization
#'
#'@description Input the reconstructed dissimilarity matrix and output the hierarchical clustering result to help judge the number of potential cell clusters.
#'
#'@param rcdm the reconstructed dissimilarity matrix
#'
#'@examples Get_hclust(rsm)
#'
#'@return  hierarchical clustering result
#'
#'@export
#'
#'@import pheatmap
Get_hclust <- function(rcdm){
  return(pheatmap(rcdm))
}

#'@title Assess clustering performace
#'
#'@description Calculate ARI to evaluate clustering performance
#'
#'@param hc_obj the hierarchical clustering object from reconstructed dissimilarity matrix
#'
#'@param clu the potential cluster number which could be obtained from hierarchical clustering heatmap of Get_hclust
#'
#'@param ref_label the available or known cell type label
#'
#'@examples
#'true_label <- load('syth_label.RData')
#'k_clu <- 6
#'ARI_test <- Comp_ARI(hc_obj,k_clu,true_label)
#
#'@return the clustering performance
#'
#'@export
#'
#'@import aricode
Comp_ARI <- function(hc_obj,clu,ref_label){
  pre    <- cutree(hc_obj,k=clu)
  return(ARI(pre,ref_label))
}
