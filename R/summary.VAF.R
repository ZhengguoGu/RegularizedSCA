#'Display a summary of the results of \code{VAF()}.
#'
#'@param object Object of class inheriting from 'VAF'.
#'@param ...  Argument to be passed to or from other methods. 
#'@examples
#'\dontrun{
#'## S3 method for class 'VAF'
#'summary(object)
#'}
#'
#'@export
summary.VAF <- function(object, ...){
  
    cat(sprintf("\nProportion of VAF for each block:\n"))
    print(object$block, digits)
    
    cat(sprintf("\nProportion of VAF for each component of each block:\n"))
    print(object$component, digits)
    
}