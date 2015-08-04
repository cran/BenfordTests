\name{usq.benftest}
\alias{usq.benftest}
\title{
Freedman-Watson U-square Test for Benford's Law
}
\description{
\code{usq.benftest} takes any numerical vector reduces the sample to the specified number of significant digits and performs the Freedman-Watson test for discreet distributions between the first digits' distribution and Benford's distribution to assert if the data conforms to Benford's law.
}
\usage{
usq.benftest(x = NULL, digits = 1, pvalmethod = "simulate", pvalsims = 10000)
}
\arguments{
  \item{x}{A numeric vector.}
  \item{digits}{An integer determining the number of first digits to use for testing, i.e. 1 for only the first, 2 for the first two etc.}
  \item{pvalmethod}{Method used for calculating the p-value. Currently only \code{"simulate"} is available.}
  \item{pvalsims}{An integer specifying the number of replicates used if \code{pvalmethod = "simulate"}.}
}
\details{
A Freedman-Watson test for discreet distributions is performed between \code{signifd(x,digits)} and \code{pbenf(digits)}.
Specifically:\deqn{U^2 = \frac{n}{9\cdot 10^{k-1}}\cdot\left[ \displaystyle\sum_{i={10^{k-1}}}^{10^{k}-2}\left( \displaystyle\sum_{j=1}^{i}(f_j^o - f_j^e) \right)^2 - \frac{1}{9\cdot 10^{k-1}}\cdot\left(\displaystyle\sum_{i={10^{k-1}}}^{10^{k}-2}\displaystyle\sum_{j=1}^{i}(f_i^o - f_i^e)\right)^2\right]
}
where \eqn{f_i^o} denotes the observed frequency of digits \eqn{i}, and \eqn{f_i^e} denotes the expected frequency of digits \eqn{i}.
\code{x} is a numeric vector of arbitrary length. Values of \code{x} should be continuous, as dictated by theory, but may also be integers.
\code{digits} should be chosen so that \code{signifd(x,digits)} is not influenced by previous rounding.
}
\value{
A list with class "\code{htest}" containing the following components:
\item{statistic }{the value of the \eqn{U^2} test statistic}
\item{p.value }{the p-value for the test}
\item{method }{a character string indicating the type of test performed}
\item{data.name	}{ a character string giving the name of the data}
}
\references{
Benford, F. (1938) The Law of Anomalous Numbers. \emph{Proceedings of the American Philosophical Society}. \bold{78}, 551--572.

Freedman, L.S. (1981) Watson's Un2 Statistic for a Discrete Distribution. \emph{Biometrika}. \bold{68}, 708--711.

Joenssen, D.W. (2013) Two Digit Testing for Benford's Law. \emph{Proceedings of the ISI World Statistics Congress, 59th Session in Hong Kong}. [available under \url{http://www.statistics.gov.hk/wsc/CPS021-P2-S.pdf}]

Watson, G.S. (1961) Goodness-of-Fit Tests on a Circle. \emph{Biometrika}. \bold{48}, 109--114.
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
#Perform Freedman-Watson U-squared Test on
#the sample's first digits using defaults
usq.benftest(X)
#p-value = 0.4847
}