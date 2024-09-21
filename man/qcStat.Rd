\name{qcStat}
\alias{qcStat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Quality control tests for eddy covariance fluxes
}
\description{
Returns the test statistics of the quality control routines described by Vitale et al (2019). 
}
\usage{
qcStat(path_rawdata, ext_tstamp=c("START", "END"), path_output=NULL, FileName=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{path_rawdata}{
  path where raw high-frequency eddy covariance data are stored.
}

  \item{ext_tstamp}{
  Character string specifying whether the timestamp in the file name is indicative of the beginning ("START") or the end ("END") of time series involved.
}

  \item{path_output}{
  path where the results will be stored. Default is NULL.
}
  \item{FileName}{
  a character string naming a file for writing. Default is NULL.
}
}
\details{
}
\value{
A data frame containing:
\item{TSTAMP}{ISO timestamp start of averaging period (format: \emph{yyyymmddHHMM}).}
\item{SADiag}{Diagnostic for the sonic anemometer.}
\item{FMR_X}{Fraction of Missing Records in raw, high-frequency, data used for X flux variable estimation.}
\item{LGD_X}{Longest Gap Duration in raw, high-frequency, data used for X flux variable estimation.}
\item{Skew_X}{Skewness of X variable.}
\item{Kurt_X}{Kurtosis of X variable.}
\item{KID0_X}{Kurtosis Index of Differenced X variable without esclusion of zero values.}
\item{KID_X}{Kurtosis Index of Differenced X variable excluding zeros, ie the effect of possible low resolution problems.}
\item{HF5_X}{Statistic of the homogeneity test applied on X variable fluctuations (percentage of data exceeding \eqn{+/- 5\sigma}).}
\item{HF10_X}{Statistic of the homogeneity test applied on X variable fluctuations (percentage of data exceeding \eqn{+/- 10\sigma}).}
\item{HD5_X}{Statistic of the homogeneity test applied on differenced X variable (percentage of data exceeding \eqn{+/- 5\sigma}).}
\item{HD10_X}{Statistic of the homogeneity test applied on differenced X variable (percentage of data exceeding \eqn{+/- 10\sigma}).}
\item{AL1_X}{Autocorrelation at lag 1.}
\item{DDI_X}{Maximum number of binned data points (number of bins proportional to the sample size).}
\item{DIP_X}{P-value of the Hartigans' Dip Test for Unimodality.}
\item{COV_XY}{Covariance between X and Y variables.}
\item{N0_X}{Number of repeated consecutive values affecting X flux variable.}
\item{LSR_X}{Statistic of the Low Signal Resolution test for X flux variable.}
\item{M98_X}{Statistic of the non-stationary ratio test by Mahrt (1988) for X flux variable.}
}
\references{
Mahrt L (1998) Flux sampling errors for aircraft and towers, J. Atmos. Ocean. Tech., 15, 416-429, \url{https://doi.org/10.1175/1520-0426(1998)015<0416:fsefaa>2.0.co;2}.

Vitale, D. Fratini, G. Bilancia, M. Nicolini, G. Sabbatini, S. Papale, D.
A robust data cleaning procedure for eddy covariance flux measurements,
Biogeosciences, 17, 1367--1391, doi: \url{https://doi.org/10.5194/bg-17-1367-2020}, 2020.
}
\author{
Domenico Vitale
}
\note{
The output of \code{\link{qcStat}} constitutes one of input files required by the \code{\link{ecworkset}} function.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
PATH_RAWDATA <- system.file("extdata", "DE-HoH_EC_201907301200_v01.csv", package = "RFlux")
QC_STAT <- qcStat(path_rawdata=PATH_RAWDATA, ext_tstamp="END", path_output=NULL, FileName=NULL)
QC_STAT
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{datagen}% __ONLY ONE__ keyword per line
