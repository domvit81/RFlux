\name{despiking}
\alias{despiking}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Despiking algorithms}
\description{
A suite of despiking algorithms for micro-meteorological variables
}
\usage{
despiking(x, mfreq, variant, wsignal, wscale, wby = 1, zth = 5, alpha = 0.01)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
	a vector of the observed time-series values.
	}
  \item{mfreq}{
  the main frequency of the observed time series (24 or 48 for hourly and half-hourly time series, respectively; 10 or 20 for raw EC data acquired at 10Hz or 20Hz, respectively).
	}
  \item{variant}{
	"v1" mainly designed for (half-)hourly EC fluxes (missing values are allowed); "v2" mainly designed for (half-)hourly meteo time series (a low percentage of missing values is allowed); "v3" mainly designed for high-frequency EC raw data.
	}
  \item{wsignal}{
  the window width used to estimate the underlying signal (only for v2 and v3 variants).
  }
  \item{wscale}{
	the window width used to estimate the local scale paramater (only for v2 and v3 variants).
	}
  \item{wby}{
  calculate the scale parameter at every by-th point rather than every point (default=1). Large values of wby reduce the computational time, but can introduce bias in the scale parameter estimation.
	}
  \item{zth}{
  the threshold value of the z-sigma rule ((only for v2 and v3 variants).
	}
  \item{alpha}{
  the significance level used in the outlier detection rule (only for v1 variant).
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
A list containing the following components:
  \item{despiked_ts }{A vector of the despiked time-series.}
  \item{spike_loc }{Integer. Location of the detected spikes.}
}
\references{
Fried, R. Schettlinger, K. Matthias Borowski, M. robfilter: Robust
  Time Series Filters. R package version  4.1.1.  \url{https://CRAN.R-project.org/package=robfilter}, 2018.

Hafen, R. stlplus: Enhanced Seasonal Decomposition of Time Series by
  Loess. R package version 0.5.1. \url{https://CRAN.R-project.org/package=stlplus}, 2016.

Vitale, D. Fratini, G. Bilancia, M. Nicolini, G. Sabbatini, S. Papale, D.
A robust data cleaning procedure for eddy covariance flux measurements,
Biogeosciences, 17, 1367--1391, doi: \url{https://doi.org/10.5194/bg-17-1367-2020}, 2020.

Vitale, D. A perfomance evaluation of despiking algorithms for eddy covariance data. Nature Scientific Reports (in review).
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
  }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods }% use one of  RShowDoc("KEYWORDS")
\keyword{robust }% use one of  RShowDoc("KEYWORDS")
