#'A heatmap of the cross-validation result
#'
#'A heatmap of mean squared prediction errors (MSPE) against Lasso and Group 
#'Lasso tuning parameters. Note that on the x-axis 
#'are Group Lasso tuning parameter values, and the on the y-axis 
#'are the Lasso tuning parameter values.
#'
#'@param x A object for plot.
#'@param ...  Argument to be passed to or from other methods. 
#'
#'@examples
#'\dontrun{
#'## S3 method for class 'CVsparseSCA'
#'plot(x)
#'}
#'
#'@export
plot.CVsparseSCA <- function(x, ...){
  
  plotlog <- x$plotlog
  

  length_lasso <- length(x$Lasso_values)
  lasso <- 1:length_lasso
  length_glasso <- length(x$Glasso_values)
  glasso <- 1:length_glasso
  
  grid <- expand.grid(ROW=lasso, COL=glasso)
  grid$HIGHT <- c(x$MSPE)
  Acol.regions <- colorspace::diverge_hsv(n=50)
  colramp <- grDevices::colorRampPalette(c('yellow',  'green',  'blue'))
  
  if(plotlog == 1){
    p <- lattice::levelplot(HIGHT~COL*ROW, grid, col.regions = colramp,
              contour=T, labels=T,
              scales = list(x=list(at = c(1, length_glasso), labels=c(round(x$Glasso_values[1], digits = 3), round(x$Glasso_values[length_glasso], digits = 3))),
                            y=list(at = c(1, length_lasso), labels=c(round(x$Lasso_values[1], digits = 3), round(log(x$Lasso_values[length_lasso]), digits = 3)))),
              xlab = "Group Lasso tuning parameters", 
              ylab = "Lasso tuning parameters (equal spaced on the log scale)",
              main = "Mean squared prediction errors")
  }else{
    p <- lattice::levelplot(HIGHT~COL*ROW, grid, col.regions = colramp,
              contour=T, labels=T,
              scales = list(x=list(at = c(1, length_glasso), labels=c(round(x$Glasso_values[1], digits = 3), round(x$Glasso_values[length_glasso], digits = 3))),
                            y=list(at = c(1, length_lasso), labels=c(round(x$Lasso_values[1], digits = 3), round(x$Lasso_values[length_lasso], digits = 3)))),
              xlab = "Group Lasso tuning parameters", 
              ylab = "Lasso tuning parameters",
              main = "Mean Squared Prediction Errors")
    
    }
  
  
return(p)
}