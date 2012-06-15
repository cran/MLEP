\name{grr}
\alias{grr}
\title{Evaluate gradient of log likelihood}
\description{
This function evaluates gradient of log likelihood of the given likelihood polynomial for the given penetrance values.
}
\usage{
grr(para, poly)
}
\arguments{
  \item{para}{
numerical vector giving penetrance parameters. They are arranged by "alpha", "beta", and "gamma".
}
  \item{poly}{
list of numerical vector giving coefficient of likelihood polynomial. Each of the vector has \code{powers} and \code{max_power} as attribute.
}
}
\value{
Vector of gradient of log likelihood of the given likelihood polynomial for the given penetrance values.
}
\references{
Sugaya Y. MLEP: Programs for exploring the maximum likelihood estimates of penetrance parameters, Submitted. 
}
\author{
Yuki Sugaya
}




