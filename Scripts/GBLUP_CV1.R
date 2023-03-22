################################################################################
#########                 Model for hybrid prediction
#########                         GBLUP 
################################################################################

## ----------------------------------------------------------------
## ---------- 1.Packages and environment
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
## ---------- 3. Building and loading the Kernels
## ----------------------------------------------------------------

##>>---------------------------- GBLUP Kernel ---------------------

#Loading matrix
load("ADMat.Rdata") #Generated through AGHMatrix package

# Parents in CA20
G_ID = as.matrix(BLUP0[,1])

#Cutting the matrix for account only genotypes with 
#Additive kernel
K_GB = A_mat[rownames(A_mat)%in%G_ID,
             colnames(A_mat)%in%G_ID]

#Dominance kernel
K_GBD =  D_mat[rownames(D_mat)%in%G_ID,
               colnames(D_mat)%in%G_ID]


#---------------- Creating the ETA for additive effect for STM/MTM (GBLUP)

A.GB=list(list(K=K_GB,model='RKHS'))

#---------------- Creating the ETA for additive + dominance effect for STM/MTM (GBLUP)

AD.GB=list(list(K=K_GB,model='RKHS'),
           list(K=K_GBD,model='RKHS'))


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
S_GB = data.frame()

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
              ETA=A.GB,
              nIter = nIter,
              burnIn = burnIn,
              thin = thin,
              saveAt='SGB_'
    )
    
    PC = cor(Y1[Pos_NA],fm$yHat[Pos_NA])
    S_GB = rbind(S_GB,data.frame(k_Fold=p,Acc=PC,Trait=colnames(dat0[j])))
    }
  }
}


## -----------------
## ---- 5.2 The multitrait model (GBLUP)
## -----------------

# Set file for accuracy storage
M_GB = data.frame()

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
                   ETA = A.GB,
                   resCov = list(type="DIAG"),
                   nIter = nIter,
                   burnIn = burnIn,
                   thin = thin,
                   saveAt='MGB_')
  
  PC = diag(cor(Y1[Pos_NA,],fm$ETAHat[Pos_NA,]))
  M_GB = rbind(M_GB,data.frame(Acc=PC))
  }
}


## ----------------------------------------------------------------
## ---------- 6. GBLUP model (A + D)
## ----------------------------------------------------------------

## -----------------
## ---- 6.1 The single trait model (GBLUP)
## -----------------


##>>----- Set file for accuracy storage
S_GBD = data.frame()

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
              ETA=AD.GB,
              nIter = nIter,
              burnIn = burnIn,
              thin = thin,
              saveAt='SGBD_'
    )
    
    PC = cor(Y1[Pos_NA],fm$yHat[Pos_NA])
    S_GBD = rbind(S_GBD,data.frame(k_Fold=p,Acc=PC,Trait=colnames(dat0[j])))
    }
  }
}


## -----------------
## ---- 6.2 The multitrait model (GBLUP)
## -----------------

# Set file for accuracy storage
M_GBD = data.frame()

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
                   ETA = AD.GB,
                   resCov = list(type="DIAG"),
                   nIter = nIter,
                   burnIn = burnIn,
                   thin = thin,
                   saveAt='MGBD_')
  
  PC = diag(cor(Y1[Pos_NA,],fm$ETAHat[Pos_NA,]))
  M_GBD = rbind(M_GBD,data.frame(Acc=PC))
  }
}


save(S_GB,S_GBD,M_GB,M_GBD, file = "Accuracy_GBLUP.Rdata")

#################################### The end #########################################


