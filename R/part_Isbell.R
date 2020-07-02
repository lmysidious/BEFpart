#' Isbell's non-spatiotemporal partition
#'
#' @param MM An vector contains the monoculture yields.
#' @param YY An vector contains the mixture yields with compatible dimensions to MM.
#' @param PP An vector contains the expected relative mixture yields with compatible dimensions to MM. The defalut is PP=array(1/dim(MM)[1],dim=dim(MM)).
#'
#' @return net biodiversity effect, complementary effect, non-random overyielding effect, selection effect
#' @export
#'
#' @examples {
#' MM = c(3,2)
#' YY = c(2.5,1)
#' part.LH(c(3,2),c(2.5,1))
#' }
#'
part.Isbell<- function(MM,YY,PP=array(1/length(MM)[1],dim=length(MM))){
  ## non-spatiotemporal version of Isbell's partition

  WW <- YY/length(YY)
  nn <- length(MM)[1]

  M <- mean(MM)
  Z <- mean(ZZ)
  X <- mean(XX)

  NBE <- sum(MM*ZZ)
  tc <-  nn*M*Z ## complementarity effect
  no <- nn*cov2(MM,ZZ-WW) ## nonrandom  effect
  se <- nn*cov2(MM,WW) ## insurance effect

  bef.partition <- c(NBE=NBE,ce=ce,no=no,se=se)
  bef.partition
}
