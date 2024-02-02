################################################################################
#########                  Model for hybrid prediction
#########                     Gaussian Kernel - CV1
################################################################################

## ----------------------------------------------------------------
## ---------- 1. Packages and env
## ----------------------------------------------------------------

rm(list=ls())

require(BGLR)
require(tidyverse)
require(MTM)
require(BMTME)

setwd("") #Set work directory


## ----------------------------------------------------------------
## ---------- 2. Organizing the pedigree information
## ----------------------------------------------------------------

##>>-----Reading the BLUP and pedigree information
BLUP0 = read.table("CA20.txt", h=TRUE) #*** For the other two environments (FL and WIS, just change the dataset)

##>>----- Splitting for later use
dat0 = BLUP0[,c(2:9)] 

colnames(dat0) = c("STAND","DTP","HP","EW","EL","TPF","PH","EH") #*** Change for other traits and environments
dat1 = as.matrix(dat0)

BLUP1 = BLUP0 %>% 
  tibble()

## ----------------------------------------------------------------
## ---------- 3. Building  the Kernels
## ----------------------------------------------------------------

# Loading SNP matrix
load("Markers_SweetHybrid.Rdata")
Markers = Markers_SweetHybrid

# For additive and dominance kernels
KGK = getKernel(Markers = Markers_SweetHybrid, 
                 method = "euclidean",
                 GenoID = as.matrix(ID_BG)
                 MM_threshold = 0.3, 
                 maf_thresh = 0.01)

## ----------------------------------------------------------------
## ---------- 3. Building  the ETA for BGLR and MTM 
## ----------------------------------------------------------------

K_GK = KGK[[1]]
K_GKD = KGK[[2]]
h = KGK[[3]]
hD = KGK[[4]]

###>>>---------------- The ETA

A.GK=list(list(K1 = exp(-h[1]*K_GK), model='RKHS'),
          list(K2 = exp(-h[2]*K_GK), model='RKHS'),
          list(K3 = exp(-h[3]*K_GK), model='RKHS'))

###>>>---- The ETA
AD.GK=list(list(K1 = exp(-h[1]*K_GK), model='RKHS'),
           list(K2 = exp(-h[2]*K_GK), model='RKHS'),
           list(K3 = exp(-h[3]*K_GK), model='RKHS'),
           list(K4 = exp(-hD[1]*K_GKD), model='RKHS'),
           list(K5 = exp(-hD[2]*K_GKD), model='RKHS'),
           list(K6 = exp(-hD[3]*K_GKD), model='RKHS'))



## ----------------------------------------------------------------
## ---------- 4. Bayesian parameters
## ----------------------------------------------------------------

nIter = 30000
burnIn = 3000
thin = 10


## ----------------------------------------------------------------
## ---------- 5. GBLUP model + A
## ----------------------------------------------------------------

## -----------------
## ---- 5.1 The single trait model (GBLUP)
## -----------------

##>>----- Set file for accuracy storage
S_GK = data.frame()

#Data for the loop
traits = c("STAND","DTP","HP","EW","EL","TPF","PH","EH") #*** Change for other traits and environments


set.seed(0928761)
for (i in 1:5){
  
##>>----- Partitions for a 5-FCV
PT_ls = k_fold(BLUP1[,1], K=5)
PT_ls = PT_ls$CrossValidation_list

##>>----- Preparing the data entry and output

for (j in 1:length(traits)){
  
  #Data for the loop
  Y1 = dat1[,j]
  
  ##>>----- Model using STM

  for(p in 1:5){
    
    Y_NA = Y1
    Pos_NA = PT_ls[[p]]
    Y_NA[Pos_NA] = NA
    
    fm = BGLR(y=Y_NA,
              ETA=A.GK,
              nIter = nIter,
              burnIn = burnIn,
              thin = thin,
              saveAt='SGB_'
    )
    
    PC = cor(Y1[Pos_NA],fm$yHat[Pos_NA])
    S_GK = rbind(S_GK,data.frame(k_Fold=p,Acc=PC,Trait=colnames(dat0[j])))
    
    }
  }
}


## -----------------
## ---- 5.2 The multitrait model (GBLUP)
## -----------------

# Set file for accuracy storage
M_GK = data.frame()

set.seed(0928761)
for (i in 1:5){
##>>----- Partitions for a 5-FCV
PT_ls = k_fold(BLUP1[,1], K=5)
PT_ls = PT_ls$CrossValidation_list

##>>----- Model using MTM
Y1 = dat1

for(p in 1:5){
  
  Y_NA = Y1
  Pos_NA = PT_ls[[p]]
  Y_NA[Pos_NA,] = NA
  
  fm <- Multitrait(y=Y_NA,
                   ETA = A.GK,
                   resCov = list(type="DIAG"),
                   nIter = nIter,
                   burnIn = burnIn,
                   thin = thin,
                   saveAt='MGK_')
  
  PC = diag(cor(Y1[Pos_NA,],fm$ETAHat[Pos_NA,]))
  M_GK = rbind(M_GK,data.frame(Acc=PC))
  }
}


## ----------------------------------------------------------------
## ---------- 6. GBLUP model (A + D)
## ----------------------------------------------------------------

## -----------------
## ---- 6.1 The single trait model (GBLUP)
## -----------------


##>>----- Set file for accuracy storage
S_GKD = data.frame()


#Data for the loop
traits = c("STAND","DTP","HP","EW","EL","TPF","PH","EH") #*** Change for other traits and environments

set.seed(0928761)
for (i in 1:5){
##>>----- Partitions for a 5-FCV
PT_ls = k_fold(BLUP1[,1], K=5)

PT_ls = PT_ls$CrossValidation_list


##>>----- Preparing the data entry and output

for (j in 1:length(traits)){
  
  #Data for the loop
  Y1 = dat1[,j]
  
  ##>>----- Model using STM

  for(p in 1:5){
    
    Y_NA = Y1
    Pos_NA = PT_ls[[p]]
    Y_NA[Pos_NA] = NA
    
    fm = BGLR(y=Y_NA,
              ETA=AD.GK,
              nIter = nIter,
              burnIn = burnIn,
              thin = thin,
              saveAt='SGBK_'
    )
    
    PC = cor(Y1[Pos_NA],fm$yHat[Pos_NA])
    S_GKD = rbind(S_GKD,data.frame(k_Fold=p,Acc=PC,Trait=colnames(dat0[j])))
    
    }
  }
}


## -----------------
## ---- 6.2 The multitrait model (GBLUP)
## -----------------

# Set file for accuracy storage
M_GKD = data.frame()

set.seed(0928761)
for (i in 1:5){
##>>----- Partitions for a 5-FCV
PT_ls = k_fold(BLUP1[,1], K=5)
PT_ls = PT_ls$CrossValidation_list

##>>----- Model using MTM
Y1 = dat1

for(p in 1:5){
  
  Y_NA = Y1
  Pos_NA = PT_ls[[p]]
  Y_NA[Pos_NA,] = NA
  
  fm <- Multitrait(y=Y_NA,
                   ETA = AD.GK,
                   resCov = list(type="DIAG"),
                   nIter = nIter,
                   burnIn = burnIn,
                   thin = thin,
                   saveAt='MGBK_')
  
  PC = diag(cor(Y1[Pos_NA,],fm$ETAHat[Pos_NA,]))
  M_GKD = rbind(M_GKD,data.frame(Acc=PC))
  } 
}


#Save an object to a file
save(S_GK, S_GKD, M_GK, M_GKD, file = "Accuracy_RKHS.Rdata")

#################################### The end #########################################


