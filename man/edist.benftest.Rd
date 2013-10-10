\name{edist.benftest}
\alias{edist.benftest}
\title{Euclidean Distance Test for Benford's Law}
\description{
\code{edist.benftest} takes any numerical vector reduces the sample to the specified number of significant digits and performs a goodness-of-fit test based on the Euclidean distance between the first digits' distribution and Benford's distribution to assert if the data conforms to Benford's law.
}
\usage{
edist.benftest(x = NULL, digits = 1, pvalmethod = "simulate", pvalsims = 10000)
}
\arguments{
  \item{x}{A numeric vector.}
  \item{digits}{An integer determining the number of first digits to use for testing, i.e. 1 for only the first, 2 for the first two etc.}
  \item{pvalmethod}{Method used for calculating the p-value. Currently only \code{"simulate"} is available.}
  \item{pvalsims}{An integer specifying the number of replicates used if \code{pvalmethod = "simulate"}.}
}
\details{
A statistical test is performed utilizing the Euclidean distance between \code{signifd(x,digits)} and \code{pbenf(digits)}.
Specifically: \deqn{d = \sqrt{\displaystyle\sum_{i=10^{k-1}}^{10^k-1}\left(f_i^o - f_i^e\right)^2}}
where \eqn{f_i^o} denotes the observed frequency of digits \eqn{i}, and \eqn{f_i^e} denotes the expected frequency of digits \eqn{i}.
\code{x} is a numeric vector of arbitrary length.
Values of \code{x} should be continuous, as dictated by theory, but may also be integers.
\code{digits} should be chosen so that \code{signifd(x,digits)} is not influenced by previous rounding.
}
\value{
A list with class "\code{htest}" containing the following components:
\item{statistic }{the value of the Euclidean distance test statistic}
\item{p.value }{the p-value for the test}
\item{method }{a character string indicating the type of test performed}
\item{data.name	}{ a character string giving the name of the data}
}
\references{
Benford, F. (1938) The Law of Anomalous Numbers. \emph{Proceedings of the American Philosophical Society}. \bold{78}, 551--572.

Cho, W.K.T. and Gaines, B.J. (2007) Breaking the (Benford) Law: Statistical Fraud Detection in Campaign Finance. \emph{The American Statistician}. \bold{61}, 218--223.

Morrow, J. (2010) \emph{Benford's Law, Families of Distributions and a Test Basis}. [available under \url{http://www.johnmorrow.info/projects/benford/benfordMain.pdf}]

}
\author{
Dieter William Joenssen \email{Dieter.Joenssen@TU-Ilmenau.de}
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
#Perform a Euclidean Distance Test on the
#sample's first digits using defaults
edist.benftest(X,pvalmethod ="simulate")
#p-value = 0.6085
}