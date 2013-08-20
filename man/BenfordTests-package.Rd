\name{BenfordTests-package}
\alias{BenfordTests-package}
\alias{BenfordTests}
\docType{package}
\title{Statistical Tests for Benford's Law.}
\description{
This package contains seven specialized statistical tests and support functions
for determining if numerical data could conform to Benford's law.
}
\details{
\tabular{ll}{
Package: \tab BenfordTests\cr
Type: \tab Package\cr
Version: \tab 1.0.2\cr
Date: \tab 2013-08-20\cr
License: \tab GPL-3\cr
}
\code{BenfordTests} is the implementation of the seven most commonly used goodness-of-fit (GOF) tests to assess if data conforms to Benford's law.
\cr
Tests include:\cr
	Pearson \eqn{\chi^2} statistic (Pearson (1900))\cr
	Kolmogorov-Smirnov \emph{D} statistic (Kolmogorov (1933))\cr
	Freedman's modification of Watson's \eqn{U^2} statistic (Freedman (1981), Watson (1961))\cr
	Chebyshev distance \emph{m} statistic (Leemis (2000))\cr
	Euclidean distance \emph{d} statistic (Cho and Gaines (2007))\cr
	Judge-Schechter mean deviation \eqn{a^*} statistic (Judge and Schechter (2009))\cr
	Joenssen's \eqn{J_P^2} statistic, a Shapiro-Francia type correlation test (Shapiro and Francia (1972))\cr
\cr
All tests may be performed using more than one leading digit.
All tests simulate the specific p-values required for statistical inference, while p-values for the \eqn{\chi^2} and \emph{D} statistics may also be determined using their asymptotic distributions.
}
\author{
Dieter William Joenssen

Maintainer: Dieter William Joenssen <Dieter.Joenssen@TU-Ilmenau.de>
}
\references{
 
 Benford, F. (1938) The Law of Anomalous Numbers. \emph{Proceedings of the American Philosophical Society}. \bold{78}, 551--572.
 
 Cho, W.K.T. and Gaines, B.J. (2007) Breaking the (Benford) Law: Statistical Fraud Detection in Campaign Finance. \emph{The American Statistician}. \bold{61}, 218--223.

 Freedman, L.S. (1981) Watson's Un2 Statistic for a Discrete Distribution. \emph{Biometrika}. \bold{68}, 708--711.
 
 Joenssen, D.W. (2013) Two Digit Testing for Benford's Law. \emph{Proceedings of the ISI World Statistics Congress, 59th Session in Hong Kong}. [available under \url{http://www.statistics.gov.hk/wsc/CPS021-P2-S.pdf}]

 Judge, G. and Schechter, L. (2009) Detecting Problems in Survey Data using Benford's Law. \emph{Journal of Human Resources}. \bold{44}, 1--24.
 
 Kolmogorov, A.N. (1933) Sulla determinazione empirica di una legge di distibuzione. \emph{Giornale dell'Istituto Italiano degli Attuari}. \bold{4}, 83--91.
 
 Leemis, L.M., Schmeiser, B.W. and Evans, D.L. (2000) Survival Distributions Satisfying Benford's law. \emph{The American Statistician}. \bold{54}, 236--241.
 
 Newcomb, S. (1881) Note on the Frequency of Use of the Different Digits in Natural Numbers. \emph{American Journal of Mathematics}. \bold{4}, 39--40.
 
 Pearson, K. (1900) On the Criterion that a Given System of Deviations from the Probable in the Case of a Correlated System of Variables is Such that it can be Reasonably Supposed to have Arisen from Random Sampling. \emph{Philosophical Magazine Series 5}. \bold{50}, 157--175.
 
 Shapiro, S.S. and Francia, R.S. (1972) An Approximate Analysis of Variance Test for Normality. \emph{Journal of the American Statistical Association}. \bold{67}, 215--216.
 
 Watson, G.S. (1961) Goodness-of-Fit Tests on a Circle. \emph{Biometrika}. \bold{48}, 109--114.

 }

\keyword{package}
\keyword{distribution}
\keyword{datagen}
\keyword{htest}
\keyword{manip}

\examples{
#Set the random seed to an arbitrary number
set.seed(421)
#Create a sample satisfying Benford's law
X<-rbenf(n=20)
#Look at sample
X
#Look at the first digits of the sample
signifd(X)

#Perform a Chi-squared Test on the sample's first digits using defaults
chisq.benftest(X)
#p-value = 0.648
}
