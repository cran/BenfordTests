\name{signifd}
\alias{signifd}
\title{
Leading Digits
}
\description{
Returns the specified number of significant digits for each element of a given vector.
}
\usage{
signifd(x = NULL, digits = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A numeric vector.}
  \item{digits}{An integer determining the number of first digits to use for testing, i.e. 1 for only the first, 2 for the first two etc.}
}
\details{
\code{x} is a numeric vector of arbitrary length. Unlike other solutions, this function will work reliably with all real numbers.
}
\value{
Returns a vector of integers the same length as the input vector \code{x}.
}

\author{
Dieter William Joenssen \email{Dieter.Joenssen@TU-Ilmenau.de}
}
\seealso{
\code{\link{chisq.benftest}}; \code{\link{ks.benftest}}; \code{\link{usq.benftest}}; \code{\link{mdist.benftest}}; \code{\link{edist.benftest}}; \code{\link{meandigit.benftest}}; \code{\link{jpsq.benftest}}
}

\keyword{ manip }

\examples{
#Set the random seed to an arbitrary number
set.seed(421)
#Create a sample satisfying Benford's law
X<-rbenf(n=20)
#Look at the first digits of the sample
signifd(X)
#should be:
#[1] 6 1 5 2 7 5 7 4 3 1 2 2 1 5 6 2 2 3 1 3
}