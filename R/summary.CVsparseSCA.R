#'Display a summary of the results of \code{cv_sparseSCA()}.
#'
#'@param object Object of class inheriting from 'CVsparseSCA'.
#'@param disp The default is \code{"tuning"}; in this case, the recommended tuning 
#'            parameter values are presented. 
#'            If \code{"estimatedPT"}, then the estimated component loading and estimated
#'            component score matrices (based on the recommended tuning paramter
#'            values) are presented.
#'            If \code{"full"}, then information is displayed regarding 1) the 
#'            recommended tuning parameter values, 2) the estimated component 
#'            loading and estimated component score matrices (based on the 
#'            recommended tuning parameter values), 3) # of variable selected,
#'            4) Mean squared prediction error (MSPE), 
#'            5) standard errors for MSPE, 6) Lasso and Group Lasso tuning 
#'            parameter values that have been evaluated.
#'@param ...  Argument to be passed to or from other methods. 
#'@examples
#'\dontrun{
#'## S3 method for class 'CVsparseSCA'
#'summary(object, disp="full")
#'}
#'
#'@export
summary.CVsparseSCA <- function(object, disp, ...){
  
  
  if(missing(disp)){
    disp <- "tuning"
  }
  
  if(disp == "tuning"){
    
    cat(sprintf("\nRecommended tuning parameter values are:\n"))
    print(object$RecommendedLambda)
    
  } else if(disp == "estimatedPT"){
    
    cat(sprintf("\nEstimated component loading matrix, given the recommended tuning parameter values are:\n"))
    print(object$P_hat)
    
    cat(sprintf("\nEstimated component score matrix, given the recommended tuning parameter values are:\n"))
    print(object$T_hat)
    
  }else if(disp == "full"){
    
    cat(sprintf("\nRecommended tuning parameter values are:\n"))
    print(object$RecommendedLambda)
    
    cat(sprintf("\nEstimated component loading matrix, given the recommended tuning parameter values are:\n"))
    print(object$P_hat)
    
    cat(sprintf("\nEstimated component score matrix, given the recommended tuning parameter values are:\n"))
    print(object$T_hat)
    
    #cat(sprintf("\nGiven each value for Group Lasso tuning parameters, the proper region for Lasso tuning parameter values are:\n"))
    #print(object$Lambdaregion)
    
    cat(sprintf("\n# of variable selected:\n"))
    print(object$VarSelected)
    
    cat(sprintf("\nMean squared prediction error (MSPE):\n"))
    print(object$MSPE)
    
    cat(sprintf("\nstandard errors for MSPE:\n"))
    print(object$SE_MSE)
    
    cat(sprintf("\nLasso tuning parameter values that have been evaluated:\n"))
    print(object$Lasso_values)
    
    cat(sprintf("\nGroup Lasso tuning parameter values that have been evaluated:\n"))
    print(object$Glasso_values)
  
  }else{
    stop("either simple or full")
  }
  
}