#' Loreau-Hector partition
#'
#' @param MM An vector contains the monoculture yields.
#' @param YY An vector contains the mixture yields with compatible dimensions to MM.
#' @param PP An vector contains the expected relative mixture yields with compatible dimensions to MM. The defalut is PP=array(1/dim(MM)[1],dim=dim(MM)).
#'
#' @return net biodiversity effect, complementary effect, selection effect
#' @export
#'
#' @examples {
#' MM = c(3,2)
#' YY = c(2.5,1)
#' part.LH(c(3,2),c(2.5,1))
#' }
part.LH <- function(MM,YY,PP=array(1/length(MM)[1],dim=length(MM))){
  ## Loreau-Hector partition
  ZZ <- YY/MM - PP ## delta-z in the text
  nn <- length(MM)[1]

  M <- mean(MM)
  Z <- mean(ZZ)

  NBE <- sum(MM*ZZ)
  ce <- nn*M*Z ## complementarity
  se <- sum((ZZ-Z)*(MM-M)) ## selection=cov2(ZZ,MM)

  bef.partition <- c(NBE=NBE,ce=ce,se=se)
  bef.partition
}
