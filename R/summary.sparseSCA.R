#'Display a summary of the results of Cross-validation of \code{sparseSCA()}.
#'
#'@param object Object of class inheriting from 'sparseSCA'.
#'@param disp The default is \code{simple}; in this case, the recommended tuning 
#'            parameter values will be displayed. 
#'            If \code{full}, then information is displayed regarding 1) the 
#'            recommended tuning parameter values, 2) the proper region 
#'            for Lasso tuning parameter values, given a Group Lasso tuning 
#'            parameter value, 3) # of variable selected,
#'            4) Predicted residual sum of squares (PRESS), 
#'            5) standard errors for PRESS, 6) Lasso and Group Lasso tuning 
#'            parameter values that have been evaluated.
#'@examples
#'\dontrun{
#'## S3 method for class 'sparseSCA'
#'summary(object, disp="full")
#'}
#'
#'@export
summary.sparseSCA <- function(obj_sSCA, disp){
  
  PRESS <- obj_sSCA$PRESS
  SE_MSE <- obj_sSCA$SE_MSE
  VarSelected <- obj_sSCA$VarSelected
  Lasso_values <- obj_sSCA$Lasso_values
  Glasso_values <- obj_sSCA$Glasso_values
  Lambdaregion <- obj_sSCA$Lambdaregion
  RecommendedLambda <- obj_sSCA$RecommendedLambda
  
  if(missing(disp)){
    disp <- "simple"
  }
  
  if(disp == "simple"){
    
    cat(sprintf("\nRecommended tuning parameter values are:\n"))
    print(RecommendedLambda)
    
  }else if(disp == "full"){
    
    cat(sprintf("\nRecommended tuning parameter values are:\n"))
    print(RecommendedLambda)
    
    cat(sprintf("\nGiven each value for Group Lasso tuning parameters, the proper region for Lasso tuning parameter values are:\n"))
    print(Lambdaregion)
    
    cat(sprintf("\n# of variable selected:\n"))
    print(VarSelected)
    
    cat(sprintf("\nPredicted residual sum of squares (PRESS):\n"))
    print(PRESS)
    
    cat(sprintf("\nstandard errors for PRESS:\n"))
    print(SE_MSE)
    
    cat(sprintf("\nLasso tuning parameter values that have been evaluated:\n"))
    print(Lasso_values)
    
    cat(sprintf("\nGroup Lasso tuning parameter values that have been evaluated:\n"))
    print(Glasso_values)
  
  }else{
    stop("either simple or full")
  }
  
}