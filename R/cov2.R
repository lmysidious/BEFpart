#' sample covariance scaled by n
#'
#' @param x a numeric vector, matrix or data frame
#' @param y a numeric vector, matrix or data frame, with compatible dimensions to x
#'
#' @return sample covariance scaled by n, rather than n-1
#' @export
#'
#' @examples cov2(c(1,2),c(2,3))
cov2 <- function(x,y){
  cov(x,y)*(length(x)-1)/length(x)
}
