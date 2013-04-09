\name{chi_square_benford}
\alias{chi_square_benford}
\title{Pearsons's Chi-squared Goodness-of-Fit Test for Benford's Law}
\description{
\code{chi_square_benford} takes any numerical vector reduces the sample to the specified number of significant digits and performs Pearson's chi-square goodness-of-fit test to assert if the data conforms to Benford's law.
}
\usage{
chi_square_benford(x = NULL, first_digits = 1, pvalmethod = "asymptotic", pvalsims = 10000)
}
\arguments{
  \item{x}{A numeric vector.}
  \item{first_digits}{An integer determining the number of first digits to use for testing, i.e. 1 for only the first, 2 for the first two etc.}
  \item{pvalmethod}{Method used for calculating the p-value. Either \code{"asymptotic"} or \code{"simulate"}.}
  \item{pvalsims}{An integer specifying the number of replicates to use if \code{pvalmethod = "simulate"}.}
}
\details{
A chi-square goodness-of fit test is performed on \code{leading_digits(x,first_digits)} versus \code{pbenf(first_digits)}. \code{x} is a numeric vector of arbitrary length. Values of \code{x} should be continuous, as dictated by theory, but may also be integers.
\code{first_digits} should be chosen so that \code{leading_digits(x,first_digits)} is not influenced by previous rounding.
}
\value{
A list with class "\code{htest}" containing the following components:
\item{statistic }{the value of the chi-squared test statistic}
\item{p.value }{the p-value for the test}
\item{method }{a character string indicating the type of test performed}
}
\references{
Benford F. The law of anomalous numbers. Proceedings of the American Philosophical Society. 1938;78:551-572.\cr

Pearson K. On the criterion that a given system of deviations from the probable in the case of a correlated system of variables is such that it can be reasonably supposed to have arisen from random sampling. Philosophical Magazine Series 5. 1900;50(302):157-175.\cr
}
\author{
Dieter William Joenssen \email{Dieter.Joenssen@TU-Ilmenau.de}
}

\seealso{
\code{\link{pbenf}}
}

\keyword{htest}
