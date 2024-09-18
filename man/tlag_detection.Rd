\name{tlag_detection}
\alias{tlag_detection}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Time-lag detection
}
\description{
Time-lag detection by means of covariance maximization (CM), pre-whitening (PW), and pre-whitening with bootstrap (PWB) procedures.}
\usage{
tlag_detection(scalar_var, 
               tsonic_var, 
               w_var,
               mfreq,
               wdt=5,
               model = "ar", 
               LAG.MAX=mfreq*10,
               lws=0,
               uws=5,
               Rboot, 
               plot.it=FALSE,
               ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{scalar_var}{
  scalar atmospheric concentration variable (e.g. CO2 concentrations).
  }
  \item{tsonic_var}{
  sonic temperature variable.
  }
  \item{w_var}{
  vertical wind speed component variable.
  }
  \item{mfreq}{
  acquisition frequency.
  }
  \item{wdt}{
  width of the temporal window for smoothing the CCF (only for PWB procedure).
  }
\item{model}{
  an ARIMA model; if provided, it is used to prewhiten both series.
  Otherwise, an AR(p) model (the p order is selected by AIC) is fitted to the x-series and used to pre-whiten
  both series.
  }
  \item{LAG.MAX}{
  maximum lag at which to calculate the ccf. Default is 10 seconds.
  }
  \item{lws}{
  (optional) lower limit of the time lag search window (in sec) for the maximum covariance approach. Default 0 sec.
  }
  \item{uws}{
  (optional) upper limit of the time lag search window (in sec) for the maximum covariance approach. Default 5 sec.
  }
  \item{Rboot}{
  number of bootstrapped samples.
  }
  \item{plot.it}{
  Logical. If TRUE, the CCFs will be displayed (default FALSE).
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A list containing the following components (to update):

\item{mcw}{Detected time lag (timesteps) by maximum covariance using vertical wind speed (MCW).}
\item{mct}{Detected time lag (timesteps) by maximum covariance using sonic temperature (MCT).}
\item{mcw_win}{Detected time lag (timesteps) by maximum covariance using vertical wind speed in a pre-specified temporal window (MCWw).}
\item{mct_win}{Detected time lag (timesteps) by maximum covariance using sonic temperature in a pre-specified temporal window (MCTw).}
\item{pww}{Detected time lag (timesteps) by pre-whitening using vertical wind speed (PWW).}
\item{pwt}{Detected time lag (timesteps) by pre-whitening using sonic temperature (PWT).} 
\item{pwb}{Detected time lag (timesteps) by pre-whitening woth bootsptrap using vertical wind speed or sonic temperature (PWB).}
\item{pwb_lci}{Lower limit of the 95pct HDI for the PWB estimate (measure of uncertainty).}
\item{pwb_uci}{Upper limit of the 95pct HDI for PWB estimate (measure of uncertainty).}
\item{cor_pww}{Correlation estimate of the peak of PWW.}
\item{cor_pwt}{Correlation estimate of the peak of PWT.}
\item{cor_pwb}{Correlation estimate of the peak of PWB.}
\item{cov_mcw}{Covariance estimate of the peak of MCW.}
\item{cov_mct}{Covariance estimate of the peak of MCT.}
\item{cov_mcw_win}{Covariance estimate of the peak of MCWw.}
\item{cov_mct_win}{Covariance estimate of the peak of MCTw.}
\item{cov_pww}{Covariance estimate of the peak of PWW.}
\item{cov_pwt}{Covariance estimate of the peak of PWT.}
\item{cov_pwb}{Covariance estimate of the peak of PWB.} 
}
\references{
Vitale, D. Fratini, G. Helfter, C. et al. A pre-whitening with block-bootstrap cross-correlation procedure for temporal alignment of data sampled by eddy covariance systems,  
Environ Ecol Stat 31, 219–-244, doi: \url{https://doi.org/10.1007/s10651-024-00615-9}, 2024
}
\author{
Domenico Vitale
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See also \code{prewhiten} function in the TSA R package.
}
\examples{
\dontrun{
library(RFlux)
library(forecast)
library(xts)
data("closed_path_rawdata")
rawdata <- closed_path_rawdata
N <- nrow(rawdata)
timestamp_orig <- strptime("0000.00", format="\%H\%M.\%OS", tz="GMT") + seq(1, N, 1);
data.xts <- xts(rawdata, order.by=timestamp_orig);

## Detecting time lag between CO2 and vertical wind speed (W) by using an AR model (default)
tlag_co2_out <- tlag_detection(scalar_var=rawdata$CO2, tsonic_var=rawdata$T_SONIC, w_var=rawdata$W,
                               model="ar", mfreq=20, wdt=5, lws=0, uws=2, Rboot=29,
                               plot.it=TRUE)
tlag_co2_out
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods }% use one of  RShowDoc("KEYWORDS")
\keyword{htest }% __ONLY ONE__ keyword per line
