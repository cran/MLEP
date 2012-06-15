\name{mlep}
\alias{mlep}
\title{Evaluate likelihood as polynomial}
\description{
This function evaluates likelihood polynomial of penetrance parameters 
for the given pedigree data.
The function returns list of 
coefficient vector of the polynomail with attributes for 
power of the parameters.
The likelihood can be computed for each pedigree 
for the following analysis to evaluate maximum likelihood estimate of penetrance parameters.
}
\usage{
mlep(pedigree, freq)
}
\arguments{
  \item{pedigree}{
data frame of pedigree data to be analyzed. \code{pedigree} has following 6 integer columns.\cr
*Pedigree ID number\cr
*ID number\cr
*Father's ID number (0 is father is not in pedigree)\cr
*Mother's ID number (0 is mother is not in pedigree)\cr
*Sex of individual: 1 = male, 2 = female\cr
*Affected status: 2 = affected, 1 = unaffected, 0 = unknown
}
  \item{freq}{
a scalar giving the disease allele frequency}
}
\details{
Likelihood of penetrance parameters for a pedigree can be explicitely described as\cr
\code{L(alpha,beta,gamma)=sum_{i,j,k}c_{i,j,k}alpha^i beta^j gamma^k},\cr
where \code{alpha, beta}, and \code{gamma} are penetrance parameters, and 
\code{c_{i,j,k}} is coefficient of the likelihood polynomial.
Powers of the penetrance parameter \code{i}, \code{j}, and \code{k} are
converted into a value, \code{i+j*max_power+k*max_power^2}, 
to save memory usage because of its uniqueness, and
assigned to the coefficients vector as \code{powers} attribute.
The \code{max_power} is also assigned as attribute and is equal to \code{N+1},
where \code{N} is the number of informative individuals for 
the affected status in the given pedigree.
This function returns list of the coefficient vector for the given whole pedigrees.
}
\value{
List of coefficient vector of likelihood polynomail.
Each of the vector has \code{powers} and \code{max_power} as attribute.
}
\references{
Sugaya Y. MLEP: Programs for exploring the maximum likelihood estimates of penetrance parameters, Submitted. 
}
\author{
Yuki Sugaya
}




\examples{
## Load example pedigree data
data(pedigree)

## Evaluate likelihood polynomial of the example pedigree data setting disease allele frequency as 0.001
mlep(pedigree, 0.001)
}

