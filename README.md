
# SCAdata

Current stage of this package

1. functions that have been finished

  1.1 TuckerCoef.R: A function for calcuating Tucker congruence coefficients (and rotations etc.)
  
  1.2 VarSelectComDistPre.R: A function for estimating P and T by means of a MM algorithm. The sparseness in P can be predifined. 
  
  1.3 VarSelectFriedman.R: A function for estimating P and T by means of the algorithm by Friedman et al. (2010) and Yuan & Lin (2006)
  
  1.4 mySTD.R: A function for standardize the given data matrix per column over the rows
  
  1.5 VAF.R: A function for calculating proportion of variance accounted for. Used for deciding the # of components in P
  
  1.6 DISCOsca.R (pstr.R, pstrLoss.R, reflexmat.R) : A function for finding out the number of common/distinctive components in P
  
2. funtions that to be finished

  2.1 a function for deciding the region for lasso tuning parameter
  
  2.2 a function for cross validation
  
3. other things planed

  3.1 test functions systematically. Although all the functions have been tested once they are finished. But they have not been tested systematically and rigoriously with multiple blocks and varying numbers of components
  
  3.2 refine the code and reducing the running time, if possible
  
  
