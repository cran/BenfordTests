\name{K_S_benford}
\alias{K_S_benford}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Kolmogorov-Smirnov Test for Benford's Law
}
\description{
\code{K_S_benford} takes any numerical vector reduces the sample to the specified number of significant digits and performs the Kolmogorov-Smirnov goodness-of-fit test to assert if the data conforms to Benford's law.
}
\usage{
K_S_benford(x = NULL, first_digits = 1, pvalmethod = "simulate", pvalsims = 10000)
}
\arguments{
  \item{x}{A numeric vector.}
  \item{first_digits}{An integer determining the number of first digits to use for testing, i.e. 1 for only the first, 2 for the first two etc.}
  \item{pvalmethod}{Method used for calculating the p-value. Currently only \code{"simulate"} is available.}
  \item{pvalsims}{An integer specifying the number of replicates used if \code{pvalmethod = "simulate"}.}
}
\details{
A Kolmogorov-Smirnov test is performed between \code{leading_digits(x,first_digits)} and\cr\code{pbenf(first_digits)}. \code{x} is a numeric vector of arbitrary length. Values of \code{x} should be continuous, as dictated by theory, but may also be integers.
\code{first_digits} should be chosen so that \code{leading_digits(x,first_digits)} is not influenced by previous rounding.
}
\value{
A list with class "\code{htest}" containing the following components:
\item{statistic }{the value of the Kolmogorov-Smirnov \emph{D} test statistic}
\item{p.value }{the p-value for the test}
\item{method }{a character string indicating the type of test performed}
}
\references{
Benford F. The law of anomalous numbers. Proceedings of the American Philosophical Society. 1938;78:551-572.\cr

Kolmogorov AN. Sulla determinazione empirica di una legge di distibuzione. Giornale dell'Istituto Italiano degli Attuari. 1933;4:83-91.\cr
}
\author{
Dieter William Joenssen \email{Dieter.Joenssen@TU-Ilmenau.de}
}
\seealso{
\code{\link{pbenf}}
}

\keyword{ htest }