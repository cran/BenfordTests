\name{Euclidean_dist_benford}
\alias{Euclidean_dist_benford}
\title{Euclidean Distance Test for Benford's Law}
\description{
\code{Euclidean_dist_benford} takes any numerical vector reduces the sample to the specified number of significant digits and performs a goodness-of-fit test based on the Euclidean distance between the first digits' distribution and Benford's distribution to assert if the data conforms to Benford's law.
}
\usage{
Euclidean_dist_benford(x = NULL, first_digits = 1, pvalmethod = "simulate", pvalsims = 10000)
}
\arguments{
  \item{x}{A numeric vector.}
  \item{first_digits}{An integer determining the number of first digits to use for testing, i.e. 1 for only the first, 2 for the first two etc.}
  \item{pvalmethod}{Method used for calculating the p-value. Currently only \code{"simulate"} is available.}
  \item{pvalsims}{An integer specifying the number of replicates used if \code{pvalmethod = "simulate"}.}
}
\details{
A statistical test is performed utilizing the Euclidean distance between \code{leading_digits(x,first_digits)} and \code{pbenf(first_digits)}. \code{x} is a numeric vector of arbitrary length. Values of \code{x} should be continuous, as dictated by theory, but may also be integers.
\code{first_digits} should be chosen so that \code{leading_digits(x,first_digits)} is not influenced by previous rounding.
}
\value{
A list with class "\code{htest}" containing the following components:
\item{statistic }{the value of the Euclidean distance test statistic}
\item{p.value }{the p-value for the test}
\item{method }{a character string indicating the type of test performed}
}
\references{
Benford F. The law of anomalous numbers. Proceedings of the American Philosophical Society. 1938;78:551-572.\cr

Cho WKT, Gaines BJ. Breaking the (Benford) law: Statistical fraud detection in campaign finance. The American Statistician. 2007;61(4):218-223.\cr

Morrow J. Benford's law, families of distributions and a test basis. 2010.\cr http://www.johnmorrow.info/projects/benford/benfordMain.pdf.
}
\author{
Dieter William Joenssen \email{Dieter.Joenssen@TU-Ilmenau.de}
}

\seealso{
\code{\link{pbenf}}
}

\keyword{ htest }
