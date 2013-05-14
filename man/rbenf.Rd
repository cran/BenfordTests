\name{rbenf}
\alias{rbenf}

\title{
Random Sample Satisfying Benford's Law
}
\description{
Returns a random sample with length \code{n} satisfying Benford's law.
}
\usage{
rbenf(n)
}
\arguments{
  \item{n}{Number of observations.
}
}

\value{
Returns a random sample with length \code{n} satisfying Benford's law.
}
\references{
Benford F. The law of anomalous numbers. Proceedings of the American Philosophical Society. 1938;78:551-572.\cr
}
\author{
Dieter William Joenssen \email{Dieter.Joenssen@TU-Ilmenau.de}
}

\seealso{
\code{\link{qbenf}}; \code{\link{pbenf}}
}
\keyword{distribution}
\keyword{datagen}

\examples{
#Set the random seed to an arbitrary number
set.seed(421)
#Create a sample satisfying Benford's law
X<-rbenf(n=20)
#Look at sample
X
#should be
# [1] 6.159420 1.396476 5.193371 2.064033 7.001284 5.006184
#7.950332 4.822725 3.386809 1.619609 2.080063 2.242473 1.944697 5.460581
#[15] 6.443031 2.662821 2.079283 3.703353 1.364175 3.354136
}