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
Version: \tab 0.5\cr
Date: \tab 2013-04-09\cr
License: \tab GPL-3\cr
}
\code{BenfordTests} is the implementation of the seven most commonly used goodness-of-fit (GOF) tests to assess if data conforms to Benford's law.
\cr
Tests include:\cr
	Pearson \emph{chi-square} statistic (Pearson (1900))\cr
	Kolmogorov-Smirnov \emph{D} statistic (Kolmogorov (1933))\cr
	Freedman's modification of Watson's \emph{U-square} statistic (Freedman (1981), Watson (1961))\cr
	Chebyshev distance \emph{m} statistic (Leemis (2000))\cr
	Euclidean distance \emph{d} statistic (Cho and Gaines (2007))\cr
	Judge-Schechter mean deviation \emph{a-star} statistic (Judge and Schechter (2009))\cr
	Joenssen's \emph{JP-square} statistic, a Shapiro-Francia type correlation test (Shapiro and Francia (1972))\cr
\cr
All tests may be performed using more than one leading digit.
All tests simulate the specific p-values required for statistical inference, while p-values for the \emph{chi-square} and \emph{D} statistics may also be determined using their asymptotic distributions.
The package version numbers are indicative of performance and breadth of functions.
}
\author{
Dieter William Joenssen

Maintainer: Dieter William Joenssen <Dieter.Joenssen@TU-Ilmenau.de>
}
\references{
 
 Benford F. The law of anomalous numbers. Proceedings of the American Philosophical Society. 1938;78:551-572.\cr
 
 Cho WKT, Gaines BJ. Breaking the (Benford) law: Statistical fraud detection in campaign finance. The American Statistician. 2007;61(4):218-223.\cr

 Freedman LS. Watson's Un2 statistic for a discrete distribution. Biometrika. 1981;68(3):708-711.\cr
 
 Judge G, Schechter L. Detecting problems in survey data using Benford's law. Journal of Human Resources. 2009;44:1-24.\cr
 
 Kolmogorov AN. Sulla determinazione empirica di una legge di distibuzione. Giornale dell'Istituto Italiano degli Attuari. 1933;4:83-91.\cr
 
 Leemis LM, Schmeiser BW, Evans DL. Survival distributions satisfying Benford's law. The American Statistician. 2000;54(4):236-241.\cr
 
 Newcomb S. Note on the frequency of use of the different digits in natural numbers. American Journal of Mathematics. 1881;4(1):39-40.\cr
 
 Pearson K. On the criterion that a given system of deviations from the probable in the case of a correlated system of variables is such that it can be reasonably supposed to have arisen from random sampling. Philosophical Magazine Series 5. 1900;50(302):157-175.\cr
 
 Shapiro SS, Francia RS. An approximate analysis of variance test for normality. Journal of the American Statistical Association. 1972;67:215-216.\cr
 
 Watson GS. Goodness-of-fit tests on a circle. Biometrika. 1961;48:109-114.\cr

 }

\keyword{distribution}
\keyword{datagen}
\keyword{htest}
\keyword{manip}

\examples{
#Create a sample satisfying Benford's law
X<-rbenf(n=20)
X
#Look at the first digits of the sample
leading_digits(X)

#Perform a Chi-squared Test on the sample's first digits 
chi_square_benford(X)
}
