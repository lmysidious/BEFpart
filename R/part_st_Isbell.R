#' Isbell's spatiotemporal partition
#'
#' @param MM An array contains the monoculture yields, organized by a 3-dimension array: 1-species, 2-time, 3-space.
#' @param YY An array contains the mixture yields organized by a 3-dimension array: 1-species, 2-time, 3-space with compatible dimensions to MM.
#' @param PP An vector contains the expected relative mixture yields with compatible dimensions to MM. The defalut is PP=array(1/dim(MM)[1],dim=dim(MM)).
#'
#' @return net biodiversity effect, total complementary, non-random overyielding effect, average selection, temporal selection, spatial selection, spatiotemporal selection.
#' @export
#'
#' @examples {
#' MM = array(c(100,100,100,100,50,50,50,50),dim=c(2,2,2))
#' YY = array(c(50,50,50,50,50,50,50,50),dim=c(2,2,2))
#' part.st.Isbell(MM,YY)
#' }
#'
part.st.Isbell <- function(MM,YY,PP=array(1/dim(MM)[1],dim=dim(MM))){
  ###### data are organized by a 3-dimension array: 1-species, 2-time, 3-space
  RR <- apply(YY,c(2,3),function(x){x/sum(x)})
  ZZ <- YY/MM - RR
  DD <- RR - PP

  Mij <- apply(MM,c(1,2),mean)
  Mik <- apply(MM,c(1,3),mean)
  Mi <- apply(MM,1,mean)
  M <- mean(MM)

  Dij <- apply(DD,c(1,2),mean)
  Dik <- apply(DD,c(1,3),mean)
  Di <- apply(DD,1,mean)
  D <- mean(DD)

  nn <- dim(MM)[1]
  tt <- dim(MM)[2]
  pp <- dim(MM)[3]

  ### Start: define some arrays for calculation
  Mij3 <- array(Mij,dim=dim(MM))
  Mik3 <- MM*NA; for(j in 1:tt){Mik3[,j,]<-Mik}
  Mi3 <- MM*NA; for(j in 1:tt){for(i in 1:pp){Mi3[,j,i]<-Mi}}
  M3 <- array(M,dim=dim(MM))

  Dij3 <- array(Dij,dim=dim(MM))
  Dik3 <- MM*NA; for(j in 1:tt){Dik3[,j,]<-Dik}
  Di3 <- MM*NA; for(j in 1:tt){for(i in 1:pp){Di3[,j,i]<-Di}}
  D3 <- array(D,dim=dim(MM))
  ### END: define some arrays for calculation

  NBE <- sum(MM*(ZZ+DD))
  tc <- nn*tt*pp*M*mean(ZZ+DD)
  as <- nn*tt*pp*cov2(Mi,Di)
  ts <- pp*sum(unlist(lapply(1:nn,function(m){cov2(Mij[m,],Dij[m,])*tt})))
  ss <- tt*sum(unlist(lapply(1:nn,function(m){cov2(Mik[m,],Dik[m,])*pp})))
  tss <- sum((MM-Mij3-Mik3+Mi3)*(DD-Dij3-Dik3+Di3))
  no <- sum((MM-M3)*(ZZ-mean(ZZ)))

  ts.partition <- c(NBE=NBE,tc=tc,no=no,as=as,ts=ts,ss=ss,tss=tss)
  ts.partition
}
