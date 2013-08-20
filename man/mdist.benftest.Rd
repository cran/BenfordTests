\name{mdist.benftest}
\alias{mdist.benftest}
\title{
Chebyshev Distance Test (maximum norm) for Benford's Law
}
\description{
\code{mdist.benftest} takes any numerical vector reduces the sample to the specified number of significant digits and performs a goodness-of-fit test based on the Chebyshev distance between the first digits' distribution and Benford's distribution to assert if the data conforms to Benford's law.
}
\usage{
mdist.benftest(x = NULL, digits = 1, pvalmethod = "simulate", pvalsims = 10000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A numeric vector.}
  \item{digits}{An integer determining the number of first digits to use for testing, i.e. 1 for only the first, 2 for the first two etc.}
  \item{pvalmethod}{Method used for calculating the p-value. Currently only \code{"simulate"} is available.}
  \item{pvalsims}{An integer specifying the number of replicates used if \code{pvalmethod = "simulate"}.}
}
\details{
A statistical test is performed utilizing the Chebyshev distance between \code{signifd(x,digits)} and \code{pbenf(digits)}.
Specifically: \deqn{m = \max\limits_{i=10^{k-1},\ldots,10^k-1}\left|f_i^o - f_i^e\right|}
where \eqn{f_i^o} denotes the observed frequency of digits \eqn{i}, and \eqn{f_i^e} denotes the expected frequency of digits \eqn{i}.
\code{x} is a numeric vector of arbitrary length.
Values of \code{x} should be continuous, as dictated by theory, but may also be integers.
\code{digits} should be chosen so that \code{signifd(x,digits)} is not influenced by previous rounding.
}
\value{
A list with class "\code{htest}" containing the following components:
\item{statistic }{the value of the Chebyshev distance (maximum norm) test statistic}
\item{p.value }{the p-value for the test}
\item{method }{a character string indicating the type of test performed}
\item{data.name	}{ a character string giving the name of the data}
}
\references{
Benford, F. (1938) The Law of Anomalous Numbers. \emph{Proceedings of the American Philosophical Society}. \bold{78}, 551--572.

Leemis, L.M., Schmeiser, B.W. and Evans, D.L. (2000) Survival Distributions Satisfying Benford's law. \emph{The American Statistician}. \bold{54}, 236--241.

Morrow, J. (2010) \emph{Benford's Law, Families of Distributions and a Test Basis}. [available under \url{http://www.johnmorrow.info/projects/benford/benfordMain.pdf}]
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
#Perform a Chebyshev Distance Test on the
#sample's first digits using defaults
mdist.benftest(X)
#p-value = 0.6421
}