\name{meandigit.benftest}
\alias{meandigit.benftest}
\title{
Judge-Schechter Mean Deviation Test for Benford's Law
}
\description{
\code{meandigit.benftest} takes any numerical vector reduces the sample to the specified number of significant digits and performs a goodness-of-fit test based on the deviation in means of the first digits' distribution and Benford's distribution to assert if the data conforms to Benford's law.
}
\usage{
meandigit.benftest(x = NULL, digits = 1, pvalmethod = "simulate", pvalsims = 10000)
}
\arguments{
  \item{x}{A numeric vector.}
  \item{digits}{An integer determining the number of first digits to use for testing, i.e. 1 for only the first, 2 for the first two etc.}
  \item{pvalmethod}{Method used for calculating the p-value. Currently only \code{"simulate"} is available.}
  \item{pvalsims}{An integer specifying the number of replicates used if \code{pvalmethod = "simulate"}.}
}
\details{
A statistical test is performed utilizing the deviation beteen the mean digit of \code{signifd(x,digits)} and \code{pbenf(digits)}. The resulting statistic is normalized to [0,1]. \code{x} is a numeric vector of arbitrary length. Values of \code{x} should be continuous, as dictated by theory, but may also be integers.
\code{digits} should be chosen so that \code{signifd(x,digits)} is not influenced by previous rounding.
}
\value{
A list with class "\code{htest}" containing the following components:
\item{statistic }{the value of the \emph{a-star} test statistic}
\item{p.value }{the p-value for the test}
\item{method }{a character string indicating the type of test performed}
\item{data.name	}{ a character string giving the name of the data}
}
\references{
Benford F. The law of anomalous numbers. Proceedings of the American Philosophical Society. 1938;78:551-572.\cr

Judge G, Schechter L. Detecting problems in survey data using Benford's law. Journal of Human Resources. 2009;44:1-24.\cr
}
\author{
Dieter William Joenssen \email{Dieter.Joenssen@TU-Ilmenau.de}
}
\seealso{
\code{\link{pbenf}}
}

\keyword{ htest }

\examples{
#Set the random seed to an arbitrary number
set.seed(421)
#Create a sample satisfying Benford's law
X<-rbenf(n=20)
#Perform a Judge-Schechter Mean Deviation Test
#on the sample's first digits using defaults
meandigit.benftest(X)
#p-value = 0.1458
}