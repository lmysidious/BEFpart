#' Isbell's spatiotemporal partition
#'
#' @param MM An array contains the monoculture yields, organized by a 3-dimension array: 1-species, 2-time/space.
#' @param YY An array contains the mixture yields organized by a 3-dimension array: 1-species, 2-time/space with compatible dimensions to MM.
#' @param PP An vector contains the expected relative mixture yields with compatible dimensions to MM. The defalut is PP=array(1/dim(MM)[1],dim=dim(MM)).
#' @param spatial An indecator of spatial/temporal partition. The defalut is spatial=F, i.e. temporal partition.
#'
#' @return net biodiversity effect, total complementary, non-random overyielding effect, average selection, temporal/spatial selection.
#' @export
#'
#' @examples {
#' MM = array(c(100,100,100,100,50,50,50,50),dim=c(2,4))
#' YY = array(c(50,50,50,50,50,50,50,50),dim=c(2,4))
#' part.x.Isbell(MM,YY)
#' }
#'
part.x.Isbell <- function(MM,YY,PP=array(1/dim(MM)[1],dim=dim(MM)),spatial=F){
  ###### data are organized by a 3-dimension array: 1-species, 2-time, 3-space
  RR <- YY/(apply(YY,2,sum))
  ZZ <- YY/MM - RR
  DD <- RR - PP

  Mi <- apply(MM,1,mean)
  M <- mean(MM)

  Di <- apply(DD,1,mean)
  D <- mean(DD)

  nn <- dim(MM)[1]
  tt <- dim(MM)[2]

  ### END: define some arrays for calculation

  NBE <- sum(MM*(ZZ+DD))
  tc <- nn*tt*M*mean(ZZ+DD)
  as <- nn*tt*cov2(Mi,Di)
  ts <- sum(unlist(lapply(1:nn,function(m){cov2(MM[m,],DD[m,])*tt})))
  no <- sum((MM-mean(MM))*(ZZ-mean(ZZ)))

  if (spatial==F){
    ts.partition <- c(NBE=NBE,tc=tc,no=no,as=as,ts=ts)
  }
  if (spatial==T){
    ts.partition <- c(NBE=NBE,tc=tc,no=no,as=as,ss=ts);
  }
  ts.partition
}
