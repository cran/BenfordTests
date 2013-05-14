\name{pbenf}
\alias{pbenf}
\title{
Distribution Function for Benford's Distribution
}
\description{
Returns the complete Benford distribution function for a given number of first digits.
}
\usage{
pbenf(digits = 1)
}
\arguments{
  \item{digits}{An integer determining the number of first digits for which the pdf is returned, i.e. 1 for 1:9, 2 for 10:99 etc.}
}

\value{
Returns an object of class "\code{table}" containing the expected density of Benford's distribution for the given number of digits.
}
\references{
Benford F. The law of anomalous numbers. Proceedings of the American Philosophical Society. 1938;78:551-572.\cr
}
\author{
Dieter William Joenssen \email{Dieter.Joenssen@TU-Ilmenau.de}
}
\seealso{
\code{\link{qbenf}}; \code{\link{rbenf}}
}

\keyword{ distribution }

\examples{
#show Benford's predictions for the frequencies of the first digit values
pbenf(1)
}