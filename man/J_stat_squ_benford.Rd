\name{J_stat_squ_benford}
\alias{J_stat_squ_benford}

\title{
Joenssen's \emph{JP-square} Test for Benford's Law
}
\description{
\code{J_stat_squ_benford} takes any numerical vector reduces the sample to the specified number of significant digits and performs a goodness-of-fit test based on the correlation between the first digits' distribution and Benford's distribution to assert if the data conforms to Benford's law.
}
\usage{
J_stat_squ_benford(x = NULL, first_digits = 1, method = "pearson", pvalmethod = "simulate", pvalsims = 10000)
}

\arguments{
  \item{x}{A numeric vector.}
  \item{first_digits}{An integer determining the number of first digits to use for testing, i.e. 1 for only the first, 2 for the first two etc.}
  \item{method}{A character string indicating which correlation coefficient is to be computed. One of "pearson" (default), "kendall", or "spearman", can be abbreviated.}
  \item{pvalmethod}{Method used for calculating the p-value. Currently only \code{"simulate"} is available.}
  \item{pvalsims}{An integer specifying the number of replicates used if \code{pvalmethod = "simulate"}.}
}
\details{
A statistical test is performed utilizing the sign-preserved squared correlation between \code{leading_digits(x,first_digits)} and \code{pbenf(first_digits)}. \code{x} is a numeric vector of arbitrary length. Values of \code{x} should be continuous, as dictated by theory, but may also be integers.
\code{first_digits} should be chosen so that \code{leading_digits(x,first_digits)} is not influenced by previous rounding.
}
\value{
A list with class "\code{htest}" containing the following components:
\item{statistic }{the value of the \emph{JP-square} test statistic}
\item{p.value }{the p-value for the test}
\item{method }{a character string indicating the type of test performed}
}
\references{
Benford F. The law of anomalous numbers. Proceedings of the American Philosophical Society. 1938;78:551-572.\cr

Joenssen DW. A new test for Benford's distribution [abstract]. In: Abstract-proceedings of the 3rd joint Statistical Meeting DAGStat, March 18-22, 2013; Freiburg, Germany.\cr

Shapiro SS, Francia RS. An approximate analysis of variance test for normality. Journal of the American Statistical Association. 1972;67:215-216.\cr

}
\author{
Dieter William Joenssen \email{Dieter.Joenssen@TU-Ilmenau.de}
}
\seealso{
\code{\link{pbenf}}
}

\keyword{ htest }