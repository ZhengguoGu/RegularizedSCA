#'Display a summary of the results of \code{VAF()}.
#'
#'@param object Object of class inheriting from 'VAF'.
#'@param digits Specify the digits.
#'@param ...  Argument to be passed to or from other methods. 
#'@examples
#'\dontrun{
#'## S3 method for class 'VAF'
#'summary(object)
#'}
#'
#'@export
summary.VAF <- function(object, digits, ...){
  
    if(is.missing(digits)){
      digits <- NULL
    }
    cat(sprintf("\nProportion of VAF for each block:\n"))
    print(object$block, digits)
    
    cat(sprintf("\nProportion of VAF for each component of each block:\n"))
    print(object$component, digits)
    
}