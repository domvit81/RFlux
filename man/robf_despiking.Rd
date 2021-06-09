\name{robf_despiking}
\alias{robf_despiking}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Despiking eddy covariance raw data}
\description{
A parallel implementation of robust functionals for despiking high-frequency eddy covariance raw data}
\usage{
parallel_despiking(x, mfreq, wsignal = mfreq*30+1, wscale=mfreq*30+1, wby = 1, zth = 5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  a vector of the observed time-series values.
}
  \item{mfreq}{
  the main frequency of the observed time series (24 or 48 for hourly and half-hourly time series, respectively; 10 or 20 for raw EC data acquired at 10Hz or 20Hz, respectively).
}
  \item{file_length}{
  the raw data file length duration in minutes (e.g. 30 or 60).
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A list containing the following components:
  \item{ts_cleaned }{A vector of the despiked time-series.}
  \item{spike_loc }{Integer. Location of the detected spikes.}
}
\references{
Fried, R. Schettlinger, K. Matthias Borowski, M. robfilter: Robust
  Time Series Filters. R package version  4.1.1.  \url{https://CRAN.R-project.org/package=robfilter}, 2018.

Vitale, D. A perfomance evaluation of despiking algorithms for eddy covariance data, 
Scientific Reports, 11, 11628, doi: \url{https://doi.org/10.1038/s41598-021-91002-y}, 2021.
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
\keyword{optimize }% use one of  RShowDoc("KEYWORDS")
