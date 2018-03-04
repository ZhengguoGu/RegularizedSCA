#'Display a summary of the results of \code{cv_structuredSCA()}.
#'
#'@param object Object of class inheriting from 'CVstructuredSCA'.
#'@param disp The default is \code{"tuning"}; in this case, the recommended tuning 
#'            parameter values for Lasso is displayed
#'            If \code{"estimatedPT"}, then the estimated component loading and
#'            component score matrices (given the recommended tuning parameter)
#'            is displayed. 
#'            If \code{"full"}, then information is displayed regarding 1) the 
#'            recommended tuning parameter values for Lasso, 2) the estimated component 
#'            loading and component score matrices, 3) the proper region 
#'            for Lasso tuning parameter values, based on the 1SE rule, 4) mean squared
#'            prediction error (MSPE), 5) Lasso tuning 
#'            parameter values that have been evaluated.
#'@param ...  Argument to be passed to or from other methods. 
#'@examples
#'\dontrun{
#'## S3 method for class 'CVstructuredSCA'
#'summary(object, disp="full")
#'}
#'
#'@export
summary.CVstructuredSCA <- function(object, disp, ...){
  
  PRESS <- object$MSPE
  LassoSequence <- object$LassoSequence
  LassoRegion <- object$LassoRegion
  RecommendedLasso <- object$RecommendedLasso 
  
  if(missing(disp)){
    disp <- "tuning"
  }
  
  if(disp == "tuning"){
    
    cat(sprintf("\nRecommended tuning parameter value for Lasso:\n"))
    print(RecommendedLasso)
    
  }else if(disp == "estimatedPT"){
    
    cat(sprintf("\nEstimated component loading matrix, given the recommended Lasso tuning parameter:\n"))
    print(object$P_hat)
    
    cat(sprintf("\nEstimated component score matrix, given the recommended Lasso tuning parameter:\n"))
    print(object$T_hat)
    
  }else if(disp == "full"){
    
    cat(sprintf("\nRecommended tuning parameter value for Lasso:\n"))
    print(RecommendedLasso)
    
    cat(sprintf("\nEstimated component loading matrix, given the recommended Lasso tuning parameter:\n"))
    print(object$P_hat)
    
    cat(sprintf("\nEstimated component score matrix, given the recommended Lasso tuning parameter:\n"))
    print(object$T_hat)
    
    cat(sprintf("\nA region for suitable Lasso tuning parameter values based on 1SE rule:\n"))
    print(LassoRegion)
    
    
    cat(sprintf("\nMean squared prediction error (MSPE):\n"))
    print(PRESS)
    
    
    cat(sprintf("\nLasso tuning parameter values that have been evaluated:\n"))
    print(LassoSequence)
    
    
  }else{
    stop("either simple or full")
  }
}