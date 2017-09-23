#'Display a summary of the results of \code{undoShrinkage()}.
#'
#'@param object Object of class inheriting from 'undoS'.
#'@param ...  Argument to be passed to or from other methods. 
#'@examples
#'\dontrun{
#'## S3 method for class 'undoS'
#'summary(object)
#'}
#'
#'@export
summary.undoS <- function(object, ...){
  
  cat(sprintf("\nThe estimated P matrix is\n"))
  print(object$Pmatrix)
  
  cat(sprintf("\nThe estimated T matrix is\n"))
  print(object$Tmatrix)
  cat(sprintf("Note: The rows of T matrix represent observational units."))
}