\name{inst_prob_test}
\alias{inst_prob_test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Detection of instrumental malfunctions, errors in eddy covariance systems
}
\description{
A set of tests aims at detecting instrumental malfunctions, errors in eddy covariance systems.
}
\usage{
inst_prob_test(x,
               detrend=FALSE,
               hz=c(10,20),
               plot=FALSE,
               var_name=c("U","V","W","T_SONIC","CO2","H2O"),
               cex.test=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  raw high frequency eddy covariance time series.
}
  \item{detrend}{
  logical. If TRUE a linear trend is removed from data (see details).
}
  \item{hz}{
  integer. The data acquisition frequency (ie 10 or 20 HZ).
}
  \item{plot}{
  logical. If TRUE a graphical representation of the results is provided. Default FALSE.
}
  \item{var_name}{
  character. One of "U","V","W","T_SONIC","CO2","H2O" variable name can be specified to enrich the graphcial output when plot=TRUE.
}
  \item{cex.test}{
  A numerical value giving the amount by which plotting test results. Default set to 1.
}
}
\details{
The detrend option should be set in agreement with the processing options adopted for fluxes estimation. In case fluxes are estimated after mean removal (i.e. block average), then detrend should be set to FALSE. In case fluxes are estimated after linear trend removal then detrend should be set to TRUE.
}
\value{
  \item{Skew}{Skewness.}
  \item{Kurt}{Kurtosis index on original data.}
  \item{KID0}{Kurtosis index on difference data.}
  \item{KID}{Kurtosis index on difference data excluding low resolution problems.}
  \item{HFx}{Homogeneity test statistic based on fluctuation.}
  \item{HDx}{Homogeneity test statistic based on difference data.}
  \item{AL1}{Autocorrelation function estimate at lag 1.}
  \item{DDI}{Data Distribution Integrity test statistic. Maximum number of data points in a bin of the histogram with bins defined according to Freedman-Diaconis rule.}
  \item{DIP}{P-value of the Hartigans' Dip test for unimodality.}
}
\references{
Vitale, D. Fratini, G. Bilancia, M. Nicolini, G. Sabbatini, S. Papale, D. A robust data cleaning procedure for eddy covariance flux measurements, Biogeosciences, 17, 1367--1391, doi: \url{https://doi.org/10.5194/bg-17-1367-2020}, 2020.

Freedman, D. and Diaconis, P. On the histogram as a density estimator: L2 theory. Zeitschrift fur Wahrscheinlichkeitstheorie und verwandte Gebiete, 57, 453--476, doi: \url{https://doi.org/10.1007/BF01025868}, 1981.

Hartigan, P.M. Computation of the Dip Statistic to Test for Unimodality. Applied Statistics (JRSS C) 34, 320--325. Corresponding Fortran code of AS 217 available from Statlib, \url{http://lib.stat.cmu.edu/apstat/217}, 1985.

Hartigan, J.A. and Hartigan, P.M. The Dip Test of Unimodality. Annals of Statistics, 13, 70--84, 1985.

}
\author{
Domenico Vitale
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
data(closed_path_rawdata)
inst_prob_test(closed_path_rawdata$W, detrend=FALSE, h=20, plot=TRUE, var_name="W")
}
\keyword{univar}% use one of  RShowDoc("KEYWORDS")
