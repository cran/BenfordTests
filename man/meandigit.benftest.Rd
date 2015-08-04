\name{meandigit.benftest}
\alias{meandigit.benftest}
\title{
Judge-Schechter Mean Deviation Test for Benford's Law
}
\description{
\code{meandigit.benftest} takes any numerical vector reduces the sample to the specified number of significant digits and performs a goodness-of-fit test based on the deviation in means of the first digits' distribution and Benford's distribution to assert if the data conforms to Benford's law.
}
\usage{
meandigit.benftest(x = NULL, digits = 1, pvalmethod = "asymptotic", pvalsims = 10000)
}
\arguments{
  \item{x}{A numeric vector.}
  \item{digits}{An integer determining the number of first digits to use for testing, i.e. 1 for only the first, 2 for the first two etc.}
  \item{pvalmethod}{Method used for calculating the p-value. Either \code{"asymptotic"} or \code{"simulate"}.}
  \item{pvalsims}{An integer specifying the number of replicates used if \code{pvalmethod = "simulate"}.}
}
\details{
A statistical test is performed utilizing the deviation between the mean digit of \code{signifd(x,digits)} and \code{pbenf(digits)}.
Specifically:\deqn{a^*=\frac{|\mu_k^o-\mu_k^e|}{\left(9\cdot10^{k-1}\right)-\mu_k^e}}
where \eqn{\mu_k^o} is the observed mean of the chosen \eqn{k} number of digits, and \eqn{\mu_k^e} is the expected/true mean value for Benford's predictions.
\eqn{a^*} conforms asymptotically to a truncated normal distribution under the null-hypothesis, i.e.,
\deqn{a^*\sim  truncnorm\left(\mu=0,\sigma=\sigma_B,a=0,b=\infty\right)}
\code{x} is a numeric vector of arbitrary length.
Values of \code{x} should be continuous, as dictated by theory, but may also be integers.
\code{digits} should be chosen so that \code{signifd(x,digits)} is not influenced by previous rounding.
}
\value{
A list with class "\code{htest}" containing the following components:
\item{statistic }{the value of the \eqn{a^*} test statistic}
\item{p.value }{the p-value for the test}
\item{method }{a character string indicating the type of test performed}
\item{data.name	}{ a character string giving the name of the data}
}
\references{
Benford, F. (1938) The Law of Anomalous Numbers. \emph{Proceedings of the American Philosophical Society}. \bold{78}, 551--572.

Judge, G. and Schechter, L. (2009) Detecting Problems in Survey Data using Benford's Law. \emph{Journal of Human Resources}. \bold{44}, 1--24.
}
\author{
Dieter William Joenssen \email{Dieter.Joenssen@googlemail.com}
}
\seealso{
\code{\link{pbenf}}, \code{\link{simulateH0}}
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