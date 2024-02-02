#########################################
#
# Package: NA
#
# File: getKernel.R
# Contains: getKernel, missing_data, maf_filter
#
# Written by Marco Antonio Peixoto
#
# First version: Mar-2022
# Last update: Dec-2024
#
# License: GPL-3
#
#########################################

#' Prediction of total genetic value for a set of crosses
#'
#' @description
#' Predicts total genetic value for a set of crosses. In this case, we followed the formulae from Falconer and Mackay (1996).
#'
#' @param Markers matrix with markers information for all candidate parents,
#' coded as 0,1,2.
#' @param GenoID genotypes IDs for filtering for environments.
#' @param estimateD  dominance Gaussian kernel should be estimated? Default is FALSE
#' @param method flag for Euclidean distance
#' @param MM_threshold markers missing threshold
#' @param maf_thresh minor allele frequency threshold
#'
#' @return RKHS Gaussian kernel
#'
#' @author Marco Antonio Peixoto, \email{marco.peixotom@@gmail.com}
#'
#' @export


getKernel = function(Markers = NULL, GenoID = NULL, estimateD = FALSE, method = "euclidean", MM_threshold = 0.9, maf_thresh = 0.01){
  
  cat('Total SNPs', ncol(Markers), '\n')
  
  MarkMis = missing_data(Markers = Markers,MM_threshold) 
  cat(ncol(Markers)-ncol(MarkMis), 'SNPs dropped due to missing data threshold of', MM_threshold, '\n' )
  
  MarkMaf = maf_filter(Markers = MarkMis,maf_thresh = maf_thresh)
  cat(ncol(MarkMis)-ncol(MarkMaf), 'SNPs dropped due to maf threshold of',maf_thresh, '\n' )
  
  cat('Total SNPs', ncol(MarkMaf), '\n')
  
  ##----- Additive
  if(method == "euclidean"){
    DistMat <- as.matrix(dist(MarkMaf, method = "euclidean"))^2
  }
  
  DistMat<-as.matrix(DistMat/mean(DistMat))
  
  if(!is.null(GenoID)){
    G_ID = as.matrix(GenoID)
    
    # Cutting the matrix for account-only genotypes with
    DistMat = DistMat[rownames(DistMat)%in%G_ID,
                 colnames(DistMat)%in%G_ID]
    
    }
    
  h=round(1/median(DistMat[row(DistMat)>col(DistMat)]),2)
  h=h*c(5,1,0.2)
    
  if(isTRUE(estimateD)){
  #------ Dominance
  MarkMafD = MarkMaf
  MarkMafD[MarkMafD == 2] <- 0
  
  if(method == "euclidean"){
    DistMatD <- as.matrix(dist(MarkMafD, method = "euclidean"))^2
  }
  
  DistMatD<-as.matrix(DistMatD/mean(DistMatD))
  
  if(!is.null(GenoID)){
    G_ID = as.matrix(GenoID)
    
    # Cutting the matrix for account-only genotypes with
    DistMatD = DistMatD[rownames(DistMatD)%in%G_ID,
                      colnames(DistMatD)%in%G_ID]
    
  }
  
  hD=round(1/median(DistMatD[row(DistMatD)>col(DistMatD)]),2)
  hD=hD*c(5,1,0.2)
  
  } else {
    MarkMafD = list()
    hD = list()
  }
  
  return(list(Kernel_Add = MarkMaf,
              Kernel_Dom = MarkMafD,
              h = h,
              hD = hD))
  
}




missing_data <- function(Markers, MM_threshold){
  
  marker.missing <- apply(Markers,2,function(x)
  {return(length(which(is.na(x)))/nrow(Markers))
  })
  
  filtered <- Markers[,which(marker.missing<MM_threshold)]
 
  MarkFinal = apply(filtered, 2, FUN=function(wna) sapply(wna, function(ina) ifelse(is.na(ina), mean(wna, na.rm=TRUE), ina)))
  
  return(MarkFinal)
}


maf_filter<-function(Markers,maf_thresh){
  freq<-colMeans(Markers, na.rm=T)/2
  maf<-freq
  maf[which(maf > 0.5)] <- 1-maf[which(maf > 0.5)]
  snps<-Markers[,which(maf>maf_thresh)]
  return(snps)
}

