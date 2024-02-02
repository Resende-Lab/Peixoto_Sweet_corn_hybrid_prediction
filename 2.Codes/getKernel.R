
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
#' @param method flag for Euclidean distance
#' @param h should h be estimated
#' @param IM_threshold individual missing threshold
#' @param MM_threshold markers missing threshold
#' @param maf_thresh minor allele frequency threshold
#'
#' @return RKHS gaussian kernel
#'
#' @author Marco Antonio Peixoto, \email{marco.peixotom@@gmail.com}
#'
#' @export


getKernel = function(Markers = NULL, method = "euclidean", h = TRUE, IM_threshold = 0.9, MM_threshold = 0.9, maf_thresh = 0.01){
  
 
Markers = missing_data(Markers = Markers,IM_threshold,MM_threshold) 

Markers = maf_filter(Markers = Markers,maf_thresh = maf_thresh)
  
if(method == "euclidean"){
  Mat_RKHS<-as.matrix(dist(Markers, method = "euclidean"))^2
  
}

Mat_RKHS<-as.matrix(Mat_RKHS/mean(Mat_RKHS))

if(h == TRUE){
h=round(1/median(Mat_RKHS[row(Mat_RKHS)>col(Mat_RKHS)]),2)
h=h*c(5,1,0.2)
  
}else{
  h = 0
}

return(list(RKHS_mat = Mat_RKHS,
            h = h))
  
}




missing_data <- function(Markers,IM_threshold,MM_threshold){
 
  # Missing
  Markers = apply(Markers, 2, FUN=function(wna) sapply(wna, function(ina) ifelse(is.na(ina), mean(wna, na.rm=TRUE), ina)))
  
  ind.missing <- apply(Markers,1,function(x){
    return(length(which(is.na(x)))/ncol(Markers))
  })
  
  marker.missing <- apply(Markers,2,function(x)
  {return(length(which(is.na(x)))/nrow(Markers))
  })
  
  filtered <- Markers[which(ind.missing<IM_threshold),which(marker.missing<MM_threshold)]
  return(filtered)
}


maf_filter<-function(Markers,maf_thresh){
  freq<-colMeans(Markers, na.rm=T)/2
  maf<-freq
  maf[which(maf > 0.5)] <- 1-maf[which(maf > 0.5)]
  snps<-Markers[,which(maf>maf_thresh)]
  return(snps)
}







