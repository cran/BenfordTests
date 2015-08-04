\name{chisq.benftest}
\alias{chisq.benftest}
\title{Pearson's Chi-squared Goodness-of-Fit Test for Benford's Law}
\description{
\code{chisq.benftest} takes any numerical vector reduces the sample to the specified number of significant digits and performs Pearson's chi-square goodness-of-fit test to assert if the data conforms to Benford's law.
}
\usage{
chisq.benftest(x = NULL, digits = 1, pvalmethod = "asymptotic", pvalsims = 10000)
}
\arguments{
  \item{x}{A numeric vector.}
  \item{digits}{An integer determining the number of first digits to use for testing, i.e. 1 for only the first, 2 for the first two etc.}
  \item{pvalmethod}{Method used for calculating the p-value. Either \code{"asymptotic"} or \code{"simulate"}.}
  \item{pvalsims}{An integer specifying the number of replicates to use if \code{pvalmethod} = "simulate".}
}
\details{
A \eqn{\chi^2} goodness-of-fit test is performed on \code{signifd(x,digits)} versus \code{pbenf(digits)}.
Specifically: \deqn{\chi^2 = n\cdot\displaystyle\sum_{i=10^{k-1}}^{10^k-1}\frac{\left(f_i^o - f_i^e\right)^2}{f_i^e}}
where \eqn{f_i^o} denotes the observed frequency of digits \eqn{i}, and \eqn{f_i^e} denotes the expected frequency of digits \eqn{i}.
\code{x} is a numeric vector of arbitrary length.
Values of \code{x} should be continuous, as dictated by theory, but may also be integers.
\code{digits} should be chosen so that \code{signifd(x,digits)} is not influenced by previous rounding.
}
\value{
A list with class "\code{htest}" containing the following components:
\item{statistic }{the value of the \eqn{\chi^2} test statistic}
\item{p.value }{the p-value for the test}
\item{method }{a character string indicating the type of test performed}
\item{data.name }{a character string giving the name of the data}
}
\references{
Benford, F. (1938) The Law of Anomalous Numbers. \emph{Proceedings of the American Philosophical Society}. \bold{78}, 551--572.

Joenssen, D.W. (2013) Two Digit Testing for Benford's Law. \emph{Proceedings of the ISI World Statistics Congress, 59th Session in Hong Kong}. [available under \url{http://www.statistics.gov.hk/wsc/CPS021-P2-S.pdf}]

Pearson, K. (1900) On the Criterion that a Given System of Deviations from the Probable in the Case of a Correlated System of Variables is Such that it can be Reasonably Supposed to have Arisen from Random Sampling. \emph{Philosophical Magazine Series 5}. \bold{50}, 157--175.
}
\author{
Dieter William Joenssen \email{Dieter.Joenssen@googlemail.com}
}

\seealso{
\code{\link{pbenf}}, \code{\link{simulateH0}}
}

\keyword{htest}

\examples{
#Set the random seed to an arbitrary number
set.seed(421)
#Create a sample satisfying Benford's law
X<-rbenf(n=20)
#Perform a Chi-squared Test on the sample's 
#first digits using defaults but determine
#the p-value by simulation
chisq.benftest(X,pvalmethod ="simulate")
#p-value = 0.6401
}