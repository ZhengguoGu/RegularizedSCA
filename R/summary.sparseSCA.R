#'Display a summary of the results of Cross-validation of \code{sparseSCA()}.
#'
#'@param object Object of class inheriting from 'sparseSCA'.
#'@param disp The default is \code{simple}; in this case, the recommended tuning parameter values will be displayed. 
#'            If \code{full}, then information is displayed regarding 1) the recommended tuning parameter values, 2) the proper region 
#'            for Lasso tuning parameter values, given a Group Lasso tuning parameter value, 3) # of variable selected,
#'            4) Predicted residual sum of squares (PRESS), 5) standard errors for PRESS, 6) Lasso and Group Lasso tuning 
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
    
    print("Recommended tuning parameter values are:")
    print(RecommendedLambda)
    
  }else if(disp == "full"){
    
    print("Recommended tuning parameter values are:")
    print(RecommendedLambda)
    
    print("Given each value for Group Lasso tuning parameters, the proper region for Lasso tuning parameter values are:")
    print(Lambdaregion)
    
    print("# of variable selected:")
    print(VarSelected)
    
    print("Predicted residual sum of squares (PRESS):")
    print(PRESS)
    
    print("standard errors for PRESS:")
    print(SE_MSE)
    
    print("Lasso tuning parameter values that have been evaluated:")
    print(Lasso_values)
    
    print("Group Lasso tuning parameter values that have been evaluated:")
    print(Glasso_values)
  
  }else{
    stop("either simple or full")
  }
  
  
  
  print()
}