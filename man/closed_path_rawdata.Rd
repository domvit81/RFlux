\name{closed_path_rawdata}
\alias{closed_path_rawdata}
\docType{data}
\title{
Raw, high-frequency, data for closed path systems
}
\description{
An example of raw, high-frequency, eddy covariance data for a closed path systems. %%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("closed_path_rawdata")}
\format{
  A data frame with 36000 observations of
  \describe{
    \item{\code{U}}{Horizontal wind component along the Gill HS-50 sonic anemometer x axis (m/s).}
    \item{\code{V}}{Horizontal wind component along the Gill HS-50 sonic anemometer y axis (m/s).}
    \item{\code{W}}{Vertical wind component along the Gill HS-50 sonic anemometer z axis (m/s).}
    \item{\code{T_SONIC}}{Temperature measured by the Gill HS-50 sonic anemometer (kelvin).}
    \item{\code{CO2}}{Carbon dioxide atmospheric concentrations (mixing ratio) measured by the LI-7200 gas analyzer (ppm).}
    \item{\code{H2O}}{Water vapor atmospheric concentrations (mixing ratio) measured by the LI-7200 gas analyzer (ppt).}
    \item{\code{SA_DIAG}}{Diagnostic flag output by the Gill HS-50 sonic anemometer (dimensionless).}
    \item{\code{GA_DIAG}}{Diagnostic flag output by the LI-7200 gas analyzer (dimensionless).}
    \item{\code{T_CELL}}{Average cell temperature of the LI-7200 gas analyzer (celsius).}
    \item{\code{T_CELL_IN}}{Temperature at the inlet of the LI-7200 gas analyzer cell (celsius).}
    \item{\code{T_CELL_OUT}}{Temperature at the outlet of the LI-7200 gas analyzer cell (celsius).}
    \item{\code{PRESS_CELL}}{Cell pressure of the LI-7200 gas analyzer (Kpa).}
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
Data are from the Integrated Carbon Observation System (ICOS) European Research Infrastructure
and accessible through the ICOS Carbon Portal \url{http://www.icos-cp.eu}.}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(closed_path_rawdata)
str(closed_path_rawdata)
}
\keyword{datasets}
