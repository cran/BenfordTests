\name{signifd.seq}
\alias{signifd.seq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Sequence of Possible Leading Digits
}
\description{
Returns a vector containing all possible significant digits for a given number of places.
}
\usage{
signifd.seq(digits = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{digits}{An integer determining the number of first digits to be returned, i.e. 1 for 1:9, 2 for 10:99 etc.}
}
\value{
Returns an integer vector.
}

\author{
Dieter William Joenssen \email{Dieter.Joenssen@TU-Ilmenau.de}
}

\examples{
signifd.seq(1)
seq(from=1,to=9)==signifd.seq(1)

signifd.seq(2)
seq(from=10,to=99)==signifd.seq(2)

}