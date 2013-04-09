\name{Chebyshev_dist_benford}
\alias{Chebyshev_dist_benford}
\title{
Chebyshev Distance Test for Benford's Law
}
\description{
\code{Chebyshev_dist_benford} takes any numerical vector reduces the sample to the specified number of significant digits and performs a goodness-of-fit test based on the Chebyshev distance between the first digits' distribution and Benford's distribution to assert if the data conforms to Benford's law.
}
\usage{
Chebyshev_dist_benford(x = NULL, first_digits = 1, pvalmethod = "simulate", pvalsims = 10000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A numeric vector.}
  \item{first_digits}{An integer determining the number of first digits to use for testing, i.e. 1 for only the first, 2 for the first two etc.}
  \item{pvalmethod}{Method used for calculating the p-value. Currently only \code{"simulate"} is available.}
  \item{pvalsims}{An integer specifying the number of replicates used if \code{pvalmethod = "simulate"}.}
}
\details{
A statistical test is performed utilizing the Chebyshev distance between \cr\code{leading_digits(x,first_digits)} and \code{pbenf(first_digits)}. \code{x} is a numeric vector of arbitrary length. Values of \code{x} should be continuous, as dictated by theory, but may also be integers.
\code{first_digits} should be chosen so that \code{leading_digits(x,first_digits)} is not influenced by previous rounding.
}
\value{
A list with class "\code{htest}" containing the following components:
\item{statistic }{the value of the Chebyshev distance test statistic}
\item{p.value }{the p-value for the test}
\item{method }{a character string indicating the type of test performed}
}
\references{
Benford F. The law of anomalous numbers. Proceedings of the American Philosophical Society. 1938;78:551-572.\cr

Leemis LM, Schmeiser BW, Evans DL. Survival distributions satisfying Benford's law. The American Statistician. 2000;54(4):236-241.\cr

Morrow J. Benford's law, families of distributions and a test basis. 2010. \cr http://www.johnmorrow.info/projects/benford/benfordMain.pdf.
}
\author{
Dieter William Joenssen \email{Dieter.Joenssen@TU-Ilmenau.de}
}

\seealso{
\code{\link{pbenf}}
}

\keyword{ htest }