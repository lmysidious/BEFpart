#' spatiotemporal Loreau-Hector partition
#'
#' @param MM An array contains the monoculture yields, organized by a 3-dimension array: 1-species, 2-time, 3-space.
#' @param YY An array contains the mixture yields organized by a 3-dimension array: 1-species, 2-time, 3-space with compatible dimensions to MM.
#' @param PP An vector contains the expected relative mixture yields with compatible dimensions to MM. The defalut is PP=array(1/dim(MM)[1],dim=dim(MM)).
#'
#' @return net biodiversity effect, total complementary, average selection, temporal selection, spatial selection, spatiotemporal selection.
#' @export
#'
#' @examples {
#' MM = array(c(100,100,100,100,50,50,50,50),dim=c(2,2,2))
#' YY = array(c(50,50,50,50,50,50,50,50),dim=c(2,2,2))
#' part.st.LH(MM,YY)
#' }
part.st.LH <- function(MM,YY,PP=array(1/dim(MM)[1],dim=dim(MM))){
  ###### data are organized by a 3-dimension array: 1-species, 2-time, 3-space
  ZZ <- YY/MM - PP

  Mij <- apply(MM,c(1,2),mean)
  Mik <- apply(MM,c(1,3),mean)
  Mi <- apply(MM,1,mean)
  M <- mean(MM)

  Zij <- apply(ZZ,c(1,2),mean)
  Zik <- apply(ZZ,c(1,3),mean)
  Zi <- apply(ZZ,1,mean)
  Z <- mean(ZZ)

  nn <- dim(MM)[1]
  tt <- dim(MM)[2]
  pp <- dim(MM)[3]

  ### Start: define some arrays for calculation
  Mij3 <- array(Mij,dim=dim(MM))
  Mik3 <- MM*NA; for(j in 1:tt){Mik3[,j,]<-Mik}
  Mi3 <- MM*NA; for(j in 1:tt){for(i in 1:pp){Mi3[,j,i]<-Mi}}
  M3 <- array(M,dim=dim(MM))

  Zij3 <- array(Zij,dim=dim(MM))
  Zik3 <- MM*NA; for(j in 1:tt){Zik3[,j,]<-Zik}
  Zi3 <- MM*NA; for(j in 1:tt){for(i in 1:pp){Zi3[,j,i]<-Zi}}
  Z3 <- array(Z,dim=dim(MM))
  ### END: define some arrays for calculation

  NBE <- sum(MM*ZZ)
  tc <- nn*tt*pp*M*Z
  as <- nn*tt*pp*cov2(Mi,Zi)
  ts <- pp*sum(unlist(lapply(1:nn,function(m){cov2(Mij[m,],Zij[m,])*tt})))
  ss <- tt*sum(unlist(lapply(1:nn,function(m){cov2(Mik[m,],Zik[m,])*pp})))
  tss <- sum((MM-Mij3-Mik3+Mi3)*(ZZ-Zij3-Zik3+Zi3))

  ts.partition <- c(NBE=NBE,tc=tc,as=as,ts=ts,ss=ss,tss=tss);
  ts.partition
}
