\name{polynomial}
\alias{polynomial}
\docType{data}
\title{
Likelihood polynomial of example pedigree data
}
\description{
Likelihood polynomial of example pedigree data.
\code{polynomial} is list of numerical vector giving 
coefficient of likelihood polynomial.
Each of the vector has \code{powers} and \code{max_power} as attribute.
}
\usage{
data(polynomial)
}

\examples{
## Load likelihood polynomial of example pedigree data
data(polynomial)

## Draw perspective plot of the likelihood polynomial
PerspPenetrance(polynomial, "gamma", 0.279, theta=-60, phi=20)

## Evaluate maximum likelihood estimate of penetrance parameters
constrOptim(c(0.9,0.8,0.1), fr, grr,
	ui=rbind(c(1,0,0),c(0,1,0),c(0,0,1),c(1,-1,0),c(0,1,-1),c(-1,0,0),c(0,-1,0),c(0,0,-1)),
	ci=c(rep(0,5),rep(-1,3)), poly=polynomial, control=list(fnscale=-1))
}
\keyword{datasets}
