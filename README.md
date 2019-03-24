# For the next update (i.e., v0.5.6) the following functions will be added:
1. variable selection using index of sparseness.
2. variable selection using stability selection


# current version:RegularizedSCA v0.5.5 - An R package for regularized simultaneous component based data integration

The package performs a true integrated analysis of multiple data blocks from multiple sources. The methods included in this package combine simultaneous component analysis methods (SCA) with regularization (such as Lasso and Group Lasso). 

To use the package, please read the vignette. An article regarding this package has been submitted. We hope that the paper will be accepted and incorporated directly in this package. 

For the latest version of 'RegularizedSCA', please go to https://github.com/ZhengguoGu/RSCA/.

#Update information
1. The function maxLGlasso() has been further improved. 
2. A few bugs in\cv_sparseSCA.R(). The bugs were associated with plotting. 
3. Add a variable selection function BIC_IS(), based on BIC and Index of Sparseness.
