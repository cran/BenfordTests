\name{pbenf}
\alias{pbenf}
\title{
Probability Mass Function for Benford's Distribution
}
\description{
Returns the complete probability mass function for Benford's distribution for a given number of first digits.
}
\usage{
pbenf(digits = 1)
}
\arguments{
  \item{digits}{An integer determining the number of first digits for which the pdf is returned, i.e. 1 for 1:9, 2 for 10:99 etc.}
}
\details{Benford's distribution has the following probability mass function:
\deqn{P(d_k)=log_{10}\left(1+ d_k^{-1} \right)}
where \eqn{d_k \in \left( 10^{k-1},10^{k-1}+1, \ldots, 10^k-1 \right) } for any chosen \eqn{k} number of digits.
}
\value{
Returns an object of class "\code{table}" containing the expected density of Benford's distribution for the given number of digits.
}
\references{
Benford, F. (1938) The Law of Anomalous Numbers. \emph{Proceedings of the American Philosophical Society}. \bold{78}, 551--572.

Joenssen, D.W. (2013) Two Digit Testing for Benford's Law. \emph{Proceedings of the ISI World Statistics Congress, 59th Session in Hong Kong}. [available under \url{http://www.statistics.gov.hk/wsc/CPS021-P2-S.pdf}]
}
\author{
Dieter William Joenssen \email{Dieter.Joenssen@googlemail.com}
}
\seealso{
\code{\link{qbenf}}; \code{\link{rbenf}}
}

\keyword{ distribution }

\examples{
#show Benford's predictions for the frequencies of the first digit values
pbenf(1)
}