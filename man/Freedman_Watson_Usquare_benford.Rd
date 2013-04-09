\name{Freedman_Watson_Usquare_benford}
\alias{Freedman_Watson_Usquare_benford}
\title{
Freedman-Watson U-squared Test for Benford's Law
}
\description{
\code{Freedman_Watson_Usquare_benford} takes any numerical vector reduces the sample to the specified number of significant digits and performs the Freedman-Watson test for discreet distributions between the first digits' distribution and Benford's distribution to assert if the data conforms to Benford's law.
}
\usage{
Freedman_Watson_Usquare_benford(x = NULL, first_digits = 1, pvalmethod = "simulate", pvalsims = 10000)
}
\arguments{
  \item{x}{A numeric vector.}
  \item{first_digits}{An integer determining the number of first digits to use for testing, i.e. 1 for only the first, 2 for the first two etc.}
  \item{pvalmethod}{Method used for calculating the p-value. Currently only \code{"simulate"} is available.}
  \item{pvalsims}{An integer specifying the number of replicates used if \code{pvalmethod = "simulate"}.}
}
\details{
A Freedman-Watson test for discreet distributions is performed between \code{leading_digits(x,first_digits)} and \code{pbenf(first_digits)}. \code{x} is a numeric vector of arbitrary length. Values of \code{x} should be continuous, as dictated by theory, but may also be integers.
\code{first_digits} should be chosen so that \code{leading_digits(x,first_digits)} is not influenced by previous rounding.
}
\value{
A list with class "\code{htest}" containing the following components:
\item{statistic }{the value of the \emph{U-square} test statistic}
\item{p.value }{the p-value for the test}
\item{method }{a character string indicating the type of test performed}
}
\references{
Benford F. The law of anomalous numbers. Proceedings of the American Philosophical Society. 1938;78:551-572.\cr

Freedman LS. Watson's Un2 statistic for a discrete distribution. Biometrika. 1981;68(3):708-711.\cr

Watson GS. Goodness-of-fit tests on a circle. Biometrika. 1961;48:109-114.\cr
}
\author{
Dieter William Joenssen \email{Dieter.Joenssen@TU-Ilmenau.de}
}
\seealso{
\code{\link{pbenf}}
}

\keyword{ htest }

