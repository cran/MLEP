\name{pedigree}
\alias{pedigree}
\docType{data}
\title{
Example pedigree data
}
\description{
Example pedigree data for the evaluation of likelihood polynomial of penetrance parameters.
}
\usage{
data(pedigree)
}
\format{
  A data frame with 32 observations on the following 6 integer variables.
  \describe{
    \item{\code{PID}}{Pedigree ID number}
    \item{\code{ID}}{ID number}
    \item{\code{FID}}{Father's ID number (0 is father is not in pedigree)}
    \item{\code{MID}}{Mother's ID number (0 is father is not in pedigree)}
    \item{\code{SEX}}{Sex of individual: 1 = male, 2 = female}
    \item{\code{AFF}}{Affected status: 2 = affected, 1 = unaffected, 0 = unknown}
  }
}

\examples{
## Load example pedigree data
data(pedigree)
}
\keyword{datasets}
