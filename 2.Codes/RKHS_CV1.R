################################################################################
#########                  Model for hybrid prediction
#########                     Gaussian Kernel
################################################################################

## ----------------------------------------------------------------
## ---------- 1.Packages and env
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
#Loading matrix
load("SNP_RKHS.Rdata") #Full SNP matrix

#---------------------------------------------------------
# SNP matrix from hybrid data
SNP_HY = hybridA_RKHS

# Parents in CA20
G_ID = as.matrix(BLUP0[,1])

#Cutting the matrix for account only genotypes with 
SNP_GKA = SNP_HY[rownames(SNP_HY)%in%G_ID,]

##>>----- DISTANCE MATRIX - Creating the kernel
DA<-as.matrix(dist(SNP_GKA, method = "euclidean"))^2
DA<-DA/mean(DA)

#Pérez and de-los-Campos (2014)
h=round(1/median(DA[row(DA)>col(DA)]),2)
h=h*c(5,1,0.2)

#------------------------------------------------------------
# SNP matrix from hybrid data
SNP_D = hybridD_RKHS

# Parents in CA20
G_ID = as.matrix(BLUP0[,1])

#Cutting the matrix for account only genotypes with 
SNP_GKD = SNP_D[rownames(SNP_D)%in%G_ID,]

##>>----- DISTANCE MATRIX - Creating the kernel (Github vignette Package: https://github.com/gdlc/BGLR-R/blob/master/inst/md/RKHS.md)
DD<-as.matrix(dist(SNP_GKD, method = "euclidean"))^2
DD<-DD/mean(DD)

#Pérez and de-los-Campos (2014)
hD=round(1/median(DD[row(DD)>col(DD)]),2)
hD=hD*c(5,1,0.2)

## ----------------------------------------------------------------
## ---------- 3. Building  the ETA for BGLR and MTM 
## ----------------------------------------------------------------


#---------------------- Creating the ETA for additive effect for STM (GK)

A.GK<-list()

for(i in 1:length(h)){
  A.GK[[i]]<-list(K=exp(-h[i]*DA),model='RKHS')
}


# #---------------------- Creating the ETA for additive + dominance effect for STM (GK)

AD.GK<-list()

AD.GK=list(list(K1 = exp(-h[1]*DA), model='RKHS'),
           list(K2 = exp(-h[2]*DA), model='RKHS'),
           list(K3 = exp(-h[3]*DA), model='RKHS'),
           list(K4 = exp(-hD[1]*DD), model='RKHS'),
           list(K5 = exp(-hD[2]*DD), model='RKHS'),
           list(K6 = exp(-hD[3]*DD), model='RKHS'))



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


