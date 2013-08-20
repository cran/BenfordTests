\name{qbenf}
\alias{qbenf}
\title{
Quantile Function for Benford's Distribution
}
\description{
Returns the complete quantile function for Benford's distribution with a given number of first digits.
}
\usage{
qbenf(digits = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{digits}{An integer determining the number of first digits for which the qdf is returned, i.e. 1 for 1:9, 2 for 10:99 etc.}
}

\value{
Returns an object of class "\code{table}" containing the expected quantile function of Benford's distribution with a given number of digits.
}
\references{
Benford, F. (1938) The Law of Anomalous Numbers. \emph{Proceedings of the American Philosophical Society}. \bold{78}, 551--572.
}
\author{
Dieter William Joenssen \email{Dieter.Joenssen@TU-Ilmenau.de}
}
\seealso{
\code{\link{pbenf}}; \code{\link{rbenf}}
}
\keyword{ distribution }

\examples{
qbenf(1)

qbenf(1)==cumsum(pbenf(1))
}