\name{mahrt}
\alias{mahrt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Stationary test for eddy covariance fluxes 
}
\description{
Performs the non-stationary ratio test described by Mahrt (1998)
}
\usage{
mahrt(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  two-column dataframe containing raw high-frequency time series of vertical wind component (i.e. W) and scalar atmospheric variable (e.g. CO2).
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns the non-stationary ratio test statistic.
}
\references{
Mahrt L (1998) Flux sampling errors for aircraft and towers, J. Atmos. Ocean. Tech., 15, 416-429, \url{https://doi.org/10.1175/1520-0426(1998)015<0416:fsefaa>2.0.co;2}.
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
PATH_RAWDATA <- system.file("extdata", "DE-HoH_EC_201907301200_v01.csv", package = "RFlux")
data(closed_path_rawdata)
mahrt(data.frame("W"=closed_path_rawdata$W,"CO2"=closed_path_rawdata$CO2))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{univar}% use one of  RShowDoc("KEYWORDS")
