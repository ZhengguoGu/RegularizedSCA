

# SCAdata

Current stage of this package

1. functions that have been finished

  1.1 TuckerCoef.R: A function for calcuating Tucker congruence coefficients (and rotations etc.) (double checked on 2016-12-30)
  
  1.2 CDfriedman.R: A function for estimating P and T by means of a MM algorithm. The sparseness in P can be predifined. 
  
  1.3 CDpre.R: A function for estimating P and T by means of the algorithm by Friedman et al. (2010) and Yuan & Lin (2006)
  
  1.4 mySTD.R: A function for standardize the given data matrix per column over the rows (double checked on 2016-12-30)
  
  1.5 VAF.R: A function for calculating proportion of variance accounted for. Used for deciding the # of components in P (double checked on 2016-12-30)
  
  1.6 DISCOsca.R (pstr.R, pstrLoss.R, reflexmat.R) : A function for finding out the number of common/distinctive components in P. pstr.R, pstrLoss.R, reflexmat.R have been double checked on 2016-12-30
  
  1.7 pca_gca.R: A function for identifying common/distinctive processes by means of PCA-GCA method. 
  
  1.8 cv_CDpre.R: A cross-validation procedure for the CDpre.R algorithm
  
  1.9 component_structure.R: A function for generating a structured matrix for the CDpre.R function. (double checked on 2016-12-30)
  
  1.10 cv_CDpreKf.R: A k-fold cross-validation procedure for the CDpre.R algorithm
  
2. funtions that to be finished

  2.1 (DONE!) a function for deciding the region for lasso tuning parameter. See crossvali.R
  
  2.2 (DONE!) a function for cross validation. See crossvali.R
  
3. other things planed

  3.0 improve the code: use Vectorize() (see http://stackoverflow.com/questions/13544594/create-a-vector-function-from-a-scalar-function for a nice example).

  3.1 test functions systematically. Although all the functions have been tested once they are finished. But they have not been tested systematically and rigoriously with multiple blocks and varying numbers of components
  
  3.2 refine the code and reducing the running time, if possible
  
  
