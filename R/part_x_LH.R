#' spatiotemporal Loreau-Hector partition
#'
#' @param MM An array contains the monoculture yields, organized by a 3-dimension array: 1-species, 2-time, 3-space.
#' @param YY An array contains the mixture yields organized by a 3-dimension array: 1-species, 2-time, 3-space with compatible dimensions to MM.
#' @param PP An vector contains the expected relative mixture yields with compatible dimensions to MM. The defalut is PP=array(1/dim(MM)[1],dim=dim(MM)).
#' @param spatial An indecator of spatial/temporal partition. The defalut is spatial=F, i.e. temporal partition.
#'
#' @return net biodiversity effect, total complementary, average selection, temporal/spatial selection.
#' @export
#'
#' @examples {
#' MM = array(c(100,100,100,100,50,50,50,50),dim=c(2,4))
#' YY = array(c(50,50,50,50,50,50,50,50),dim=c(2,4))
#' part.x.LH(MM,YY)
#' }
part.x.LH <- function(MM,YY,PP=array(1/dim(MM)[1],dim=dim(MM)),spatial=F){
  ###### data are organized by a 2-dimension array: 1-species, 2-time/space
  ZZ <- YY/MM - PP

  Mi <- apply(MM,1,mean)
  M <- mean(MM)

  Zi <- apply(ZZ,1,mean)
  Z <- mean(ZZ)

  nn <- dim(MM)[1]
  tt <- dim(MM)[2]

  ### END: define some arrays for calculation

  NBE <- sum(MM*ZZ)
  tc <- nn*tt*M*Z
  as <- nn*tt*cov2(Mi,Zi)
  ts <- sum(unlist(lapply(1:nn,function(m){cov2(MM[m,],ZZ[m,])*tt})))
  if (spatial==F){
    ts.partition <- c(NBE=NBE,tc=tc,as=as,ts=ts)
  }
  if (spatial==T){
    ts.partition <- c(NBE=NBE,tc=tc,as=as,ss=ts);
  }
  ts.partition
}
