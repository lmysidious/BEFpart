% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/part_Isbell.R
\name{part.Isbell}
\alias{part.Isbell}
\title{Isbell's non-spatiotemporal partition}
\usage{
part.Isbell(MM, YY, PP = array(1/length(MM)[1], dim = length(MM)))
}
\arguments{
\item{MM}{An vector contains the monoculture yields.}

\item{YY}{An vector contains the mixture yields with compatible dimensions to MM.}

\item{PP}{An vector contains the expected relative mixture yields with compatible dimensions to MM. The defalut is PP=array(1/dim(MM)[1],dim=dim(MM)).}
}
\value{
net biodiversity effect, complementary effect, non-random overyielding effect, selection effect
}
\description{
Isbell's non-spatiotemporal partition
}
\examples{
{
MM = c(3,2)
YY = c(2.5,1)
part.LH(c(3,2),c(2.5,1))
}

}
