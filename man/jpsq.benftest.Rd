\name{jpsq.benftest}
\alias{jpsq.benftest}

\title{
Joenssen's \emph{JP-square} Test for Benford's Law
}
\description{
\code{jpsq.benftest} takes any numerical vector reduces the sample to the specified number of significant digits and performs a goodness-of-fit test based on the correlation between the first digits' distribution and Benford's distribution to assert if the data conforms to Benford's law.
}
\usage{
jpsq.benftest(x = NULL, digits = 1, pvalmethod = "simulate", pvalsims = 10000)
}

\arguments{
  \item{x}{A numeric vector.}
  \item{digits}{An integer determining the number of first digits to use for testing, i.e. 1 for only the first, 2 for the first two etc.}
  \item{pvalmethod}{Method used for calculating the p-value. Currently only \code{"simulate"} is available.}
  \item{pvalsims}{An integer specifying the number of replicates used if \code{pvalmethod = "simulate"}.}
}
\details{
A statistical test is performed utilizing the sign-preserved squared correlation between \cr\code{signifd(x,digits)} and \code{pbenf(digits)}.
Specifically: \deqn{J_P^2=sgn\left(cor\left(f^o, f^e\right)\right)\cdot cor\left(f^o, f^e\right) ^2}
where \eqn{f^o} denotes the observed frequencies and \eqn{f^e} denotes the expected frequency of digits \cr \eqn{10^{k-1},10^{k-1}+1,\ldots,10^k-1}.
\code{x} is a numeric vector of arbitrary length.
Values of \code{x} should be continuous, as dictated by theory, but may also be integers.
\code{digits} should be chosen so that \code{signifd(x,digits)} is not influenced by previous rounding.
}
\value{
A list with class "\code{htest}" containing the following components:
\item{statistic }{the value of the \eqn{J_P^2} test statistic}
\item{p.value }{the p-value for the test}
\item{method }{a character string indicating the type of test performed}
\item{data.name	}{ a character string giving the name of the data}
}
\references{
Benford, F. (1938) The Law of Anomalous Numbers. \emph{Proceedings of the American Philosophical Society}. \bold{78}, 551--572.

Joenssen, D.W. (2013) A New Test for Benford's Distribution. In: \emph{Abstract-Proceedings of the 3rd Joint Statistical Meeting DAGStat, March 18-22, 2013}; Freiburg, Germany.

Joenssen, D.W. (2013) Two Digit Testing for Benford's Law. \emph{Proceedings of the ISI World Statistics Congress, 59th Session in Hong Kong}. [available under \url{http://www.statistics.gov.hk/wsc/CPS021-P2-S.pdf}]

Shapiro, S.S. and Francia, R.S. (1972) An Approximate Analysis of Variance Test for Normality. \emph{Journal of the American Statistical Association}. \bold{67}, 215--216.

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
#Perform Joenssen's \emph{JP-square} Test
#on the sample's first digits using defaults
jpsq.benftest(X)
#p-value = 0.3241
}