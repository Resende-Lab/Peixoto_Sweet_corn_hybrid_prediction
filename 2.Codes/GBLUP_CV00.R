################################################################################
#########                 Model for CV00 and A+D
#########                         GBLUP 
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

##>>-----Reading the BLUP: This dataset contains only the individuals not shared in between 2020 and 2021
##> You can build it up based on the two datasets that are into the main folder
BLUP20 = read.table("../1.phen_files/California_20.txt", h=TRUE) #*** For the other two environments (FL and WIS, just change the dataset)

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

##>>---------------------------- Creating A+D matrices ---------------------
# Loading SNP matrix

load("Markers_SweetHybrid") # load the SNP matrix
Markers = Markers_SweetHybrid

# Generate relationship matrix using AGHMatrix package
# Additive
G_mat <- Gmatrix(Markers,
                 missingValue = "NA",
                 integer = FALSE,
                 thresh.missing = .3,
                 maf = 0.01)

# Dominance
D_mat <- Gmatrix(Markers,
                 missingValue = "NA",
                 integer = FALSE,
                 thresh.missing = .3,
                 maf = 0.01,
                 method="Vitezica")

# Parents in BG20
ID_CA20 = as.matrix(BLUP20[,1])

# Parents in BG21
ID_CA21 = as.matrix(BLUP21[,1])

# Combining BG20+BG21 and unique
ID_CA = rbind(ID_CA20,ID_CA21)

ID_CA = unique(ID_CA)

#Cutting the matrix for account-only genotypes with
#Additive kernel
K_GB = A_mat[rownames(D_mat)%in%ID_CA,
             colnames(D_mat)%in%ID_CA]

#Dominance kernel
K_GBD =  D_mat[rownames(D_mat)%in%ID_CA,
               colnames(D_mat)%in%ID_CA]

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
##>>>---- Creating design matrix for lines
#matrix
Z_L=model.matrix(~0+IDGen,data=dat0)

#Expanding G for lines additive
K_G=Z_L%*%K_GB%*%t(Z_L)

#Expanding G for lines
K_GD=Z_L%*%K_GBD%*%t(Z_L)

###>>>---------------- The ETA
A.GB=list(list(model='RKHS',K=K_G))

###>>>---------------- The ETA
AD.GB=list(add = list(model='RKHS',K=K_G),
           dom = list(model='RKHS',K=K_GD))

## ----------------------------------------------------------------
## ---------- 7. Single trait model with A
## ----------------------------------------------------------------
###>>>---- Model

#accuracy storage
S_GB = data.frame()

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
              ETA=A.GB,
              nIter = nIter, 
              burnIn = burnIn,
              thin = thin)
  
  PC = cor(Y1[Pos_NA],mod1$yHat[Pos_NA])
  S_GB = rbind(S_GB,data.frame(Acc=PC,Trait=colnames(dat0[i+1])))
  
}

## ----------------------------------------------------------------
## ---------- 8. The multitrait model (GBLUP) with A
## ----------------------------------------------------------------
###>>>---- Model

#accuracy storage
M_GB = data.frame()

#Each trait at once
Y1=as.matrix(dat0[,2:9])

#Data from 2021 to NA 
Y_NA = Y1
Pos_NA = p
Y_NA[Pos_NA,] = NA

set.seed(0928761)
#model
mod2 = Multitrait(y = Y_NA,
                  ETA=A.GB,
                  resCov = list(type="DIAG"),
                  nIter = nIter,
                  burnIn = burnIn,
                  thin = thin,
                  saveAt = 'MGB_')


PC = diag(cor(Y1[Pos_NA,],mod2$ETAHat[Pos_NA,]))
M_GB = rbind(M_GB,data.frame(Acc=PC))


## ----------------------------------------------------------------
## ---------- 7. GBLUP single trait model (A + D)
## ----------------------------------------------------------------
###>>>---- Model

#accuracy storage
S_GBD = data.frame()

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
              ETA=AD.GB,
              nIter = nIter, 
              burnIn = burnIn,
              thin = thin)
  
  PC = cor(Y1[Pos_NA],mod3$yHat[Pos_NA])
  S_GBD = rbind(S_GBD,data.frame(Acc=PC,Trait=colnames(dat0[i+1])))
  
}

## ----------------------------------------------------------------
## ---------- 8. The multitrait model (GBLUP) with A + D
## ----------------------------------------------------------------
###>>>---- Model

#accuracy storage
M_GBD = data.frame()

#Each trait at once
Y1=as.matrix(dat0[,2:9])


#Data from 2021 to NA 
Y_NA = Y1
Pos_NA = p
Y_NA[Pos_NA,] = NA

set.seed(0928761)
#model
mod4 = Multitrait(y = Y_NA,
                  ETA=AD.GB,
                  resCov = list(type="DIAG"),
                  nIter = nIter,
                  burnIn = burnIn,
                  thin = thin,
                  saveAt = 'MGB_')


PC = diag(cor(Y1[Pos_NA,],mod4$ETAHat[Pos_NA,]))
M_GBD = rbind(M_GBD,data.frame(Acc=PC))


#------------Saving the outputs
save(S_GB,S_GBD,M_GB,M_GBD, file = "Accuracy_GBLUP.Rdata")

#################################### The end #########################################


