\name{signifd.analysis}
\alias{signifd.analysis}
\title{
Graphical Analysis of First Significant Digits
}
\description{
\code{signifd.analysis} takes any numerical vector reduces the sample to the specified number of significant digits. The (relative) frequencies are then plotted so that a subjective analysis may be performed.
}
\usage{
signifd.analysis(x = NULL, digits = 1, graphical_analysis = TRUE, freq = FALSE, 
alphas = 20, tick_col = "red", ci_col = "darkgreen", ci_lines = c(.05))
}
\arguments{
  \item{x}{A numeric vector.}
  \item{digits}{An integer determining the number of first digits to use for testing, i.e. 1 for only the first, 2 for the first two etc.}
  \item{graphical_analysis}{Boolean value indicating if results should be plotted.}
  \item{freq}{Boolean value indicating if absolute frequencies should be used.}
  \item{alphas}{Either a vector containing the significance levels([0,1]) that will be shaded, or an integer defining the number of evenly spaced confidence intervals.}
  \item{tick_col}{Color code or name that will be passed to "\code{points}" for plotting.}
  \item{ci_col}{Color code or name that will be passed to "\code{polygon}" for shading the different confidence intervals. May be more than one color.}
  \item{ci_lines}{Boolean or fractional value(s) indicating significance levels where lines are drawn}
  
}
\details{
Confidence intervals are calculated from the normal distribution with \eqn{\mu_i = np_i} and \eqn{\sigma^2 = np_i(1-p_i)}, where i represents the considered digit. Be aware that the normal approximation only holds for "large" n.
}
\value{
A list containing the following components:
\item{summary }{the summary printed below the graph, a matrix of digits, their (relative) frequencies and individual p-values}
\item{CIs }{confidence intervals used for plotting as defined by parameter "\code{ci_lines}" or "\code{alphas}" if \code{ci_lines==FALSE}}
}
\references{
Benford, F. (1938) The Law of Anomalous Numbers. \emph{Proceedings of the American Philosophical Society}. \bold{78}, 551--572.

Freedman, L.S. (1981) Watson's Un2 Statistic for a Discrete Distribution. \emph{Biometrika}. \bold{68}, 708--711.
}
\author{
Dieter William Joenssen \email{Dieter.Joenssen@googlemail.com}
}
\seealso{
\code{\link{pbenf}}
}

\keyword{ hplot }

\examples{
#Set the random seed to an arbitrary number
set.seed(421)
#Create a sample satisfying Benford's law
X<-rbenf(n=20)
#Analyze the first digits using the the defaults
signifd.analysis(X)
#Turn off plot
signifd.analysis(X,graphical_analysis=FALSE)
#Use absolute frequencies
signifd.analysis(X,graphical_analysis=FALSE,freq=TRUE)
#Use five evenly spaced confidence intervals, no lines
#alphas is used for shadeing
signifd.analysis(X,graphical_analysis=TRUE,alphas=5,freq=TRUE,ci_lines=FALSE)
#Use fifty evenly spaced, gray confidence intervals, blue ticks, and lines at 
#the 1 and 5 percent confidence intervals
signifd.analysis(X,graphical_analysis=TRUE,alphas=50,freq=TRUE,tick_col="blue",
ci_col="gray",ci_lines=c(.01,.05))
}