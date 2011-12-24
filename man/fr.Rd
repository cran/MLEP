\name{fr}
\alias{fr}
\title{Evaluate log likelihood value}
\description{
This function evaluates log likelihood value of the given likelihood polynomial for the given penetrance values.
}
\usage{
fr(para, poly, minimize = F)
}
\arguments{
  \item{para}{
numerical vector giving penetrance parameters. They are arranged by "alpha", "beta", and "gamma".
}
  \item{poly}{
list of numerical vector giving coefficient of likelihood polynomial. Each of the vector has \code{powers} and \code{max_power} as attribute.
}
  \item{minimize}{
if \code{TRUE}, minus log likelihood value is returned. This option may be used to minimize the likelihood polynomial.
}
}

\value{
Scalar giving log likelihood value of the given likelihood polynomial for the given penetrance values.
}
\references{
Sugaya Y. MLEP: Programs for exploring the maximum likelihood estimates of penetrance parameters, In preparation. 
}
\author{
Yuki Sugaya
}




\examples{
## Load likelihood polynomial of example pedigree data
data(polynomial)

## Evaluate maximum likelihood estiamte of penetrance parameters by using constrOptim function
constrOptim(c(0.9,0.8,0.1), fr, grr,
	ui=rbind(c(1,0,0),c(0,1,0),c(0,0,1),c(1,-1,0),c(0,1,-1),c(-1,0,0),c(0,-1,0),c(0,0,-1)),
	ci=c(rep(0,5),rep(-1,3)),
	poly=polynomial, minimize=TRUE)
}

