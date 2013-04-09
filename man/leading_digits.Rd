\name{leading_digits}
\alias{leading_digits}
\title{
Leading Digits
}
\description{
Returns the specified number of significant digits for each element of a given vector.
}
\usage{
leading_digits(x = NULL, number = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A numeric vector.}
  \item{number}{An integer determining the number of first digits to use for testing, i.e. 1 for only the first, 2 for the first two etc.}
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
\code{\link{chi_square_benford}}; \code{\link{K_S_benford}}; \code{\link{Freedman_Watson_Usquare_benford}}; \code{\link{Chebyshev_dist_benford}}; \code{\link{Euclidean_dist_benford}}; \code{\link{J_S_avg_dev_benford}}; \code{\link{J_stat_squ_benford}}
}

\keyword{ manip }