\name{tlag_detection}
\alias{tlag_detection}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Time-lag detection
}
\description{
Optimal time-lag detection by means of the pre-whitening procedure.}
\usage{
tlag_detection(scalar_var, 
               tsonic_var, 
               w_var,
               mfreq, 
               model = "ar", 
               LAG.MAX=mfreq*10,
               lws=0,
               ups=5,
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
  \item{model}{
  an ARIMA model; if provided, it is used to prewhiten both series.
  Otherwise, an AR(p) model (the p order is selected by AIC) is fitted to the x-series and used to pre-whiten
  both series.
  }
  \item{LAG.MAX}{
  maximum lag at which to calculate the ccf. Default is 10 seconds.
  }
  \item{lws}{
  (optional) lower limit of the time lag search window for the maximum covariance approach.
  }
  \item{uws}{
  (optional) upper limit of the time lag search window for the maximum covariance approach.
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
  \item{ccf}{Output from the ccf function on the prewhitened data.}
  \item{ar}{The AR model fit to the x-series, or x.model if it is provided.}
  \item{opt_tlag }{Optimal time lag detected (timestep).}
  \item{tlag_pw }{Time lag of the main peak of the cross-correlation function after pre-whitening.}
  \item{tlag_lmax }{Local maxima of the cross-covariance function closest to tlag_pw.}
  \item{tlag_lmin }{Local minima of the cross-covariance function closest to tlag_pw.}
  \item{corr_est }{Estimated correlation at .}
  \item{cv1pct }{Critical (abs) value of the ccf at 0.01 level.}
  \item{cv5pct }{Critical (abs) value of the ccf at 0.05 level.}
}
\references{
Vitale et al. (in prep) Optimal Time Lag Detection for Eddy Covariance Data Acquisition Systems%% ~put references to the literature/web site here ~
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

## Detecting time lag betweem CO2 and vertical wind speed (W) by using an AR model (default)
tlag_co2_out <- tlag_detection(scalar_var=rawdata$CO2, tsonic_var=rawdata$T_SONIC, w_var=rawdata$W,
                               model="ar", method="boot", mfreq=20, Rboot=29, cM=.5, j=100,
                               plot.it=TRUE)
tlag_co2_out$tlag_ref
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods }% use one of  RShowDoc("KEYWORDS")
\keyword{htest }% __ONLY ONE__ keyword per line
