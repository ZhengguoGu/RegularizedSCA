#'Display a summary of the results of \code{BIC_IS()}.
#'
#'@param object Object of class inheriting from 'BIC_IS'.
#'@param disp The default is \code{"tuning"}; in this case, the recommended tuning 
#'            parameter values based on Index of Sparseness are presented. 
#'            If \code{"tuning_BIC"}, then the recommended tuning 
#'            parameter values based on BIC are presented
#'            If \code{"full"}, then recommended tuning parameter values based on Index 
#'            of Sparseness and BIC are presented. 
#'@param ...  Argument to be passed to or from other methods. 
#'@examples
#'\dontrun{
#'## S3 method for class 'BIC_IS'
#'summary(object, disp="tuning")
#'}
#'
#'@export
summary.BIC_IS <- function(object, disp, ...){
  
  
  if(missing(disp)){
    disp <- "tuning"
  }
  
  if(disp == "tuning"){
    
    cat(sprintf("\nRecommended tuning parameter values based on Index of Sparseness are:\n"))
    print(object$IS_tuning)
    
  } else if(disp == "tuning_BIC"){
    
    cat(sprintf("\nRecommended tuning parameter values based on BIC are:\n"))
    print(object$BIC_tuning)
    
    
  }else if(disp == "full"){
    
    cat(sprintf("\nRecommended tuning parameter values based on Index of Sparseness are:\n"))
    print(object$IS_tuning)
    
    cat(sprintf("\nRecommended tuning parameter values based on BIC are:\n"))
    print(object$BIC_tuning)
    
  }
  
}