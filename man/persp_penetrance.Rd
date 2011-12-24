\name{persp_penetrance}
\alias{persp_penetrance}
\title{
Draw perspective plot of log likelihood surface 
}
\description{
This function draws a perspective plot of log likelihood surface 
of the given likelihood polynomial.
The picture is drawn by 
fixing one of the given parameters.
}
\usage{
persp_penetrance(poly, fixed, value, ...)
}
\arguments{
  \item{poly}{
list of numerical vector giving coefficient of likelihood polynomail. Each of the vector has \code{powers} and \code{max_power} as attribute. Details are described in \code{mlep}.} 
  \item{fixed}{
character giving which parameter is fixed. \code{fixed} accepts only "\code{alpha}", "\code{beta}", or "\code{gamma}".
}
  \item{value}{
a positive scalar giving the \code{fixed} parameter value.
}
  \item{\dots}{
additional arguments to \code{persp}.
}
}


\references{
Sugaya Y. MLEP: Programs for exploring the maximum likelihood estimates of penetrance parameters, In preparation. 
}
\author{
Yuki Sugaya
}

\seealso{
\code{persp}
}
\examples{
## Load likelihood polynomial of example pedigree data 
data(polynomial)

## Draw perspective plot of the polynimial
persp_penetrance(polynomial, "gamma", 0.279, theta=-60, phi=20)

}

