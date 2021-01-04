\name{tlag_detection}
\alias{tlag_detection}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Time-lag detection
}
\description{
Optimal time-lag detection by means of the pre-whitening procedure.}
\usage{
tlag_detection(x, y, mfreq, x.model = ar.res, 
               AIC = FALSE, LAG.MAX = mfreq*10,
               show.plot = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  first component series (e.g. CO2 concentrations).
  }
  \item{y}{
  second component series (e.g. vertical wind speed).
  }
  \item{mfreq}{
  acquisition frequency.
  }
  \item{x.model}{
  an ARIMA model; if provided, it is used to prewhiten both series.
  Otherwise, an AR model is fitted to the x-series and used to pre-whiten
  both series.
  }
  \item{AIC}{
  Logical. If an AR model is used for prewhitening, the AR order is chosen
  by minimizing the Akaike Information Criteria (AIC=TRUE), otherwise the
  AR model is chosen as 10*log10(N) where N is the number of observations.
  Default is FALSE.
  }
  \item{LAG.MAX}{
  maximum lag at which to calculate the ccf. Default is 10 seconds.
  }
  \item{show.plot}{
  Logical. If TRUE, the CCFs will be displayed (default FALSE).
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A list containing the following components:
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

## Detecting time lag betweem H2O and vertical wind speed (W) by using an AR model (default)
tlag_h2o_out <- tlag_detection(x=rawdata$H2O, y=rawdata$W, mfreq=20, 
                               show.plot=TRUE)
tlag_h2o_out$opt_tlag
tlag_h2o_out$corr_est

## Detecting time lag between CO2 and Sonic Temperature by using an ARIMA model
mod <- auto.arima(rawdata$CO2, stationary=FALSE)
tlag_co2_out <- tlag_detection(x=rawdata$CO2, y=rawdata$T_SONIC, mfreq=20, 
                               show.plot=TRUE, x.model=mod)
tlag_co2_out$opt_tlag
tlag_co2_out$corr_est

## Time series alignment
tlag_co2 <- tlag_co2_out$opt_tlag
tlag_h2o <- tlag_h2o_out$opt_tlag
timestamp_co2 <- strptime("0000.00", format="\%H\%M.\%OS", tz="GMT") +
                 seq(-tlag_co2+1, N - tlag_co2, 1);
timestamp_h2o <- strptime("0000.00", format="\%H\%M.\%OS", tz="GMT") +
                 seq(-tlag_h2o+1, N - tlag_h2o, 1);

co2s <- xts(rawdata$CO2, order.by=timestamp_co2);
h2os <- xts(rawdata$H2O, order.by=timestamp_h2o);
rawdata_aligned <- na.omit(data.frame(merge(data.xts[,c(1:4)], 
                  co2s, h2os, data.xts[,c(7:12)], tzone="GMT")));
colnames(rawdata_aligned) <- c("U", "V", "W", "T_SONIC", "CO2", "H2O", 
                               "SA_DIAG", "GA_DIAG", 
                               "T_CELL", "T_CELL_IN", "T_CELL_OUT", "PRESS_CELL")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods }% use one of  RShowDoc("KEYWORDS")
\keyword{htest }% __ONLY ONE__ keyword per line
