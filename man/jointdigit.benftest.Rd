\name{jointdigit.benftest}
\alias{jointdigit.benftest}

\title{
A Hotelling \emph{T-square} Type Test for Benford's Law
}
\description{
\code{jointdigit.benftest} takes any numerical vector reduces the sample to the specified number of significant digits and performs a Hotelling T-square type goodness-of-fit test to assert if the data conforms to Benford's law.
}
\usage{
jointdigit.benftest(x = NULL, digits = 1, eigenvalues="all", tol = 1e-15, 
					pvalmethod = "asymptotic", pvalsims = 10000)
}

\arguments{
  \item{x}{A numeric vector.}
  \item{digits}{An integer determining the number of first digits to use for testing, i.e. 1 for only the first, 2 for the first two etc.}
  \item{eigenvalues}{How are the eigenvalues, which are used in testing, selected.}
  \item{tol}{Tolerance in detecting values that are essentially zero.}
  \item{pvalmethod}{Method used for calculating the p-value. Currently only \code{"asymptotic"} is available.}
  \item{pvalsims}{An integer specifying the number of replicates used if \code{pvalmethod = "simulate"}.}
}
\details{
A Hotelling \eqn{T^2} type goodness-of-fit test is performed on \code{signifd(x,digits)} versus \code{pbenf(digits)}.
\code{x} is a numeric vector of arbitrary length.
\strong{argument:} \code{eigenvalues} can be defined as:
\itemize{
  \item \emph{numeric}, a vector containing which eigenvalues should be used
  \item \emph{string length = 1}, eigenvalue selection scheme:
  \itemize{
  \item \emph{"all"}, use all non-zero eigenvalues
  \item \emph{"kaiser"}, use all eigenvalues larger than the mean of all non-zero eigenvalues
  }
}
Values of \code{x} should be continuous, as dictated by theory, but may also be integers.
\code{digits} should be chosen so that \code{signifd(x,digits)} is not influenced by previous rounding.
}
\value{
A list with class "\code{htest}" containing the following components:
\item{statistic }{the value of the \eqn{T^2} test statistic}
\item{p.value }{the p-value for the test}
\item{method }{a character string indicating the type of test performed}
\item{data.name	}{ a character string giving the name of the data}
\item{eigenvalues_tested   }{ a vector containing the index numbers of the eigenvalues used in testing.}
\item{eigen_val_vect   }{ the eigen values and vectors of the null distribution. computed using \code{eigen}.}
}
\references{
Benford, F. (1938) The Law of Anomalous Numbers. \emph{Proceedings of the American Philosophical Society}. \bold{78}, 551--572.

Hotelling, H. (1931). The generalization of Student's ratio. \emph{Annals of Mathematical Statistics}. \bold{2}, 360--378.

}
\author{
Dieter William Joenssen \email{Dieter.Joenssen@googlemail.com}
}
\seealso{
\code{\link{pbenf}}
}

\keyword{ htest }

\examples{
#Set the random seed to an arbitrary number
set.seed(421)
#Create a sample satisfying Benford's law
X<-rbenf(n=20)
#Perform  Test
#on the sample's first digits using defaults
jointdigit.benftest(X)
#p-value = 0.648
#Perform  Test
#using only the two largest eigenvalues
jointdigit.benftest(x=X,eigenvalues=1:2)
#p-value = 0.5176
#Perform  Test
#using the kaiser selection criterion
jointdigit.benftest(x=X,eigenvalues="kaiser")
#p-value = 0.682
}