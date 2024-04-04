# Peixoto_Sweet_Corn_hybrid_prediction

This is a repository with the scripts of the paper: Utilizing Genomic Prediction to Boost Hybrid Performance in a Sweet Corn Breeding Program (Peixoto et al.).

### Specifications

The two folders are related to the pipeline of hybrid prediction for a sweet corn dataset. It combines the phenotypic data, genotypic data, and scripts for running genomic models and cross-validations (CV). We describe each one in the next section.

*** 

### Contents

**1.Datasets**: Contains the phenotypic and the genotypic datasets for running the analyses.

1. CA20: data on BLUEs and BLUPs for the traits from the site California/2020.
2. CA21: data on BLUEs and BLUPs for the traits from the site California/2021.
3. FL20: data on BLUEs and BLUPs for the traits from the site Florida/2020.
4. FL20: data on BLUEs and BLUPs for the traits from the site Florida/2021.
5. WI20: data on BLUEs and BLUPs for the traits from the site Wisconsin/2020.
6. WI20: data on BLUEs and BLUPs for the traits from the site Wisconsin/2021.  
7. Markers_SweetHybrid: Markers coded 0,1,2 for the hybrids assessed in the locations above.  
  
**2.Scripts**: This folder has the scripts for all three different cross-validation schemes used:  

i. **CV1**: untested hybrids in tested environments. CV1 was used in all sites from 2020 (CA20, FL20, WI20).  
ii. **CV0**: tested hybrids in untested environments.  CV0 was implemented in all traits in common between 2020 and 2021 at the same site (i.e., CA20 -> CA21)  
ii. **CV00**: untested hybrids in untested environments.  CV0 was implemented in all traits in common between 2020 and 2021 at the same site (i.e., CA20 -> CA21).  
iv. **Across Env**: Implementation of Jarquin et al. (2014), section 2.5.3 from the paper. Only the traits EL, EH, and TPF were assessed. CV0 and CV00 were implemented in this case.  
v: **getKernel()** function. Used to build the distance matrix among genotypes based on Euclidian distance. Also, filter by missing values and MAF. The output can be used to build the Gaussian kernel (for more details, please look at the main scripts.

***

Any questions about the analyses, please, contact me!

Marco


Marco Antonio Peixoto  
Email: deamorimpeixotom@ufl.edu  
Page: https://marcopxt.github.io/  


