\name{simulateH0}
\alias{simulateH0}

\title{
Function for Simulating the H0-Distributions needed for BenfordTests
}
\description{
\code{simulateH0} is a wrapper function that calculates the specified test statistic under the null hypothesis a certain number of times.
}
\usage{
simulateH0(teststatistic="chisq", n=10, digits=1, pvalsims=10)
}

\arguments{
  \item{teststatistic}{Which test statistic should be used: "chisq", "edist", "jpsq", "ks", "mdist", "meandigit", or "usq".}
  \item{n}{Sample size of interest.}
  \item{digits}{An integer determining the number of first digits to use for testing, i.e. 1 for only the first, 2 for the first two etc.}
  \item{pvalsims}{An integer specifying the number of replicates to be used in simulation.}
}
\details{
Wrapper function that directly outputs the distributions of the specified test statistic under the null hypothesis. 
}
\value{
A vector of length equal to "\code{pvalsims}".
}
\references{
Benford, F. (1938) The Law of Anomalous Numbers. \emph{Proceedings of the American Philosophical Society}. \bold{78}, 551--572.

Joenssen, D.W. (2013) Two Digit Testing for Benford's Law. \emph{Proceedings of the ISI World Statistics Congress, 59th Session in Hong Kong}. [available under \url{http://www.statistics.gov.hk/wsc/CPS021-P2-S.pdf}]
}
\author{
Dieter William Joenssen \email{Dieter.Joenssen@TU-Ilmenau.de}
}
\seealso{
\code{\link{pbenf}}, \code{\link{chisq.benftest}}, \code{\link{edist.benftest}}, \code{\link{jpsq.benftest}}, \code{\link{ks.benftest}}, \code{\link{mdist.benftest}}, \\ \code{\link{meandigit.benftest}}, \code{\link{usq.benftest}}
}

\keyword{ htest }
\keyword{ distribution }

\examples{
#Set the random seed to an arbitrary number
set.seed(421)

#calculate critical value for chisquare test via simulation
quantile(simulateH0(teststatistic="chisq", n=100,digits=1,pvalsims=100000),probs=.95)

#calculate the "real" critical value
qchisq(.95,df=8)

#alternatively look at critical values for the jpsq statistic
#for different sample sizes (notice the low value for pvalsims)
set.seed(421)
apply(sapply((1:9)*10,FUN=simulateH0,teststatistic="jpsq", digits=1, pvalsims=100),
MARGIN=2,FUN=quantile,probs=.05)
}