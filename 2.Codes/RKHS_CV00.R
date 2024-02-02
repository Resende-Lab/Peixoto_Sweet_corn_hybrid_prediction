################################################################################
#########                 Model for CA21 and A+D
#########                         GK 
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
##------------------2020---------------------------

##>>-----Reading the BLUP: This dataset contains only the individuals not shared between 2020 and 2021
##> You can build it up based on the two datasets that are in the main folder
BLUP20 = read.table("California_20.txt", h=TRUE) #*** For the other two environments (FL and WIS, change the dataset)

dat20 = BLUP20

##------------------2021---------------------------
##>>-----Reading the BLUP and pedigree information
BLUP21 = read.table("CA21.txt", h=TRUE) #*** For the other two environments (FL and WIS, just change the dataset)

##>>----- Splitting for later use

dat21 = BLUP21[,c(1:9)] 

colnames(dat21) = c("IDGen","STAND","DTP","HP","EW","EL","TPF","PH","EH") #*** Change for other traits and environments

##------------------ Combining 2020-2021 ---------------------------
dat0 = rbind(dat20,dat21)

## ----------------------------------------------------------------
## ---------- 3. Building and loading the Kernels
## ----------------------------------------------------------------

# Hybrids in BG20
ID_BG20 = as.matrix(BLUP20[,1])

# Hybrids in BG21
ID_BG21 = as.matrix(BLUP21[,1])

# Combining BG20+BG21 and unique
ID_BG = rbind(ID_BG20,ID_BG21)

ID_BG = unique(ID_BG)

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
## ---------- 4. Bayesian parameters
## ----------------------------------------------------------------
nIter = 30000
burnIn = 3000
thin = 10

#NA for 2021
p <-list()
p = c(212:(211+39))


## ----------------------------------------------------------------
## ---------- 5. Creating the incidence matrix and ETA - additive
## ----------------------------------------------------------------
K_GK = KGK[[1]]
K_GKD = KGK[[2]]
h = KGK[[3]]
hD = KGK[[4]]

##>>>---- Creating design matrix for lines
Z_L=model.matrix(~0+IDGen,data=dat0)

#Expanding G for lines additive
K_G=Z_L%*%K_GK%*%t(Z_L)

#Expanding G for lines
K_GD=Z_L%*%K_GKD%*%t(Z_L)

###>>>---------------- The ETA

A.GK=list(list(K1 = exp(-h[1]*K_G), model='RKHS'),
          list(K2 = exp(-h[2]*K_G), model='RKHS'),
          list(K3 = exp(-h[3]*K_G), model='RKHS'))

###>>>---- The ETA
AD.GK=list(list(K1 = exp(-h[1]*K_G), model='RKHS'),
           list(K2 = exp(-h[2]*K_G), model='RKHS'),
           list(K3 = exp(-h[3]*K_G), model='RKHS'),
           list(K4 = exp(-hD[1]*K_GD), model='RKHS'),
           list(K5 = exp(-hD[2]*K_GD), model='RKHS'),
           list(K6 = exp(-hD[3]*K_GD), model='RKHS'))


## ----------------------------------------------------------------
## ---------- 6. Single trait model with A
## ----------------------------------------------------------------
###>>>---- Model

#accuracy storage
S_GK = data.frame()

#Traits
traits = c("STAND","DTP","HP","EW","EL","TPF","PH","EH") #*** Change for other traits and environments

set.seed(0928761)
###>>>----- Loop through the traits
for (i in 1:length(traits)){
  
  #Each trait at once
  Y1=as.matrix(dat0[,(1+i)])
  
  
  #Data from 2021 to NA 
  Y_NA = Y1
  Pos_NA = p
  Y_NA[Pos_NA] = NA
  
  
  
  #model
  mod1 = BGLR(y = Y_NA, 
              ETA=A.GK,
              nIter = nIter, 
              burnIn = burnIn,
              thin = thin)
  
  PC = cor(Y1[Pos_NA],mod1$yHat[Pos_NA])
  S_GK = rbind(S_GK,data.frame(Acc=PC,Trait=colnames(dat0[i+1])))
  
}

## ----------------------------------------------------------------
## ---------- 7. The multitrait model (GBLUP) with A
## ----------------------------------------------------------------
###>>>---- Model

#accuracy storage
M_GK = data.frame()

#Each trait at once
Y1=as.matrix(dat0[,2:9])

#Data from 2021 to NA 
Y_NA = Y1
Pos_NA = p
Y_NA[Pos_NA,] = NA

set.seed(0928761)
#model
mod2 = Multitrait(y = Y_NA,
                  ETA=A.GK,
                  resCov = list(type="DIAG"),
                  nIter = nIter,
                  burnIn = burnIn,
                  thin = thin,
                  saveAt = 'MGB_')


PC = diag(cor(Y1[Pos_NA,],mod2$ETAHat[Pos_NA,]))
M_GK = rbind(M_GK,data.frame(Acc=PC))


## ----------------------------------------------------------------
## ---------- 8. GBLUP single trait model (A + D)
## ----------------------------------------------------------------
###>>>---- Model

#accuracy storage
S_GKD = data.frame()

#Traits
traits = c("STAND","DTP","HP","EW","EL","TPF","PH","EH") #*** Change for other traits and environments

set.seed(0928761)
###>>>----- Loop through the traits
for (i in 1:length(traits)){
  
  #Each trait at once
  Y1=as.matrix(dat0[,(1+i)])
  
  
  #Data from 2021 to NA 
  Y_NA = Y1
  Pos_NA = p
  Y_NA[Pos_NA] = NA
  
  #model
  mod3 = BGLR(y = Y_NA, 
              ETA=AD.GK,
              nIter = nIter, 
              burnIn = burnIn,
              thin = thin)
  
  PC = cor(Y1[Pos_NA],mod3$yHat[Pos_NA])
  S_GKD = rbind(S_GKD,data.frame(Acc=PC,Trait=colnames(dat0[i+1])))
  
}

## ----------------------------------------------------------------
## ---------- 9. The multitrait model (GBLUP) with A + D
## ----------------------------------------------------------------
###>>>---- Model

#accuracy storage
M_GKD = data.frame()

#Each trait at once
Y1=as.matrix(dat0[,2:9])


#Data from 2021 to NA 
Y_NA = Y1
Pos_NA = p
Y_NA[Pos_NA,] = NA


set.seed(0928761)
#model
mod4 = Multitrait(y = Y_NA,
                  ETA=AD.GK,
                  resCov = list(type="DIAG"),
                  nIter = nIter,
                  burnIn = burnIn,
                  thin = thin,
                  saveAt = 'MGB_')


PC = diag(cor(Y1[Pos_NA,],mod4$ETAHat[Pos_NA,]))
M_GKD = rbind(M_GKD,data.frame(Acc=PC))


#-----------Saving the accuracy
save(S_GK,S_GKD,M_GK,M_GKD, file = "Accuracy_GK.Rdata")

#################################### The end #########################################


