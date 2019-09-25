\name{ecmd_table}
\alias{ecmd_table}
\docType{data}
\title{
Eddy covariance metadata table
}
\description{
A collection of metadata useful for the proper settings of eddy covariance data processing options.
}
\usage{data("ecmd_table")}
\format{
  A data frame with the following variables:
  \describe{
    \item{\code{DATE_OF_VARIATION_DB}}{Optional. ISO timestamp of the variation in central database. Required format: \emph{yyyymmddHHMM}.}
    \item{\code{DATE_OF_VARIATION_EF}}{ISO timestamp of effective date of variation. Required format: \emph{yyyymmdd} or \emph{yyyymmddHHMM}.}
    \item{\code{SITEID}}{Character string indicative of the EC Site's ID (\emph{CC-Xxx}).}
    \item{\code{LATITUDE}}{Latitude in decimal degree}
    \item{\code{LONGITUDE}}{Longitude in decimal degree}
    \item{\code{ALTITUDE}}{Altitude in meters}
    \item{\code{CANOPY_HEIGHT}}{Canopy height in meters}
    \item{\code{SA_MANUFACTURER}}{Sonic anemometer manufacturer (see EddyPro manual.)}
    \item{\code{SA_MODEL}}{Sonic anemometer model (see EddyPro manual.)}
    \item{\code{SA_SW_VERSION}}{Sonic anemometer software version.}
    \item{\code{SA_WIND_DATA_FORMAT}}{Format of the wind components (see EddyPro manual.)}
    \item{\code{SA_NORTH_ALIGNEMENT}}{Sonic anemometer north alignment (see EddyPro manual.)}
    \item{\code{SA_HEIGHT}}{Sonic anemometer height (see EddyPro manual.)}
    \item{\code{SA_NORTH_OFFSET}}{Sonic anemometer north offset (see EddyPro manual.)}
    \item{\code{SA_NORTH_MAGDEC}}{Sonic anemometer magnetic declination (see EddyPro manual.)}
    \item{\code{SA_INVALID_WIND_SECTOR_c1}}{Center of the first invalid wind sector.}
    \item{\code{SA_INVALID_WIND_SECTOR_w1}}{Width of the first invalid wind sector.}
    \item{\code{SA_INVALID_WIND_SECTOR_c2}}{Center of the second invalid wind sector.}
    \item{\code{SA_INVALID_WIND_SECTOR_w2}}{Width of the second invalid wind sector.}
    \item{\code{SA_INVALID_WIND_SECTOR_c3}}{Center of the third invalid wind sector.}
    \item{\code{SA_INVALID_WIND_SECTOR_w3}}{Width of the third invalid wind sector.}
    \item{\code{GA_PATH}}{Specify the eddy covariance path system: closed or open.}
    \item{\code{GA_MANUFACTURER}}{Gas analyzer manufacturer (see EddyPro manual).}
    \item{\code{GA_MODEL}}{Gas analyzer model (see EddyPro manual).}
    \item{\code{GA_SW_VERSION}}{Gas analyzer software version (see EddyPro manual).}
    \item{\code{GA_NORTHWARD_SEPARATION}}{The distance between the center of the sample volume (or the inlet of the intake tube) of the current gas analyzer and the sonic anemometer, as measured horizontally along the north-south axis (see EddyPro manual).)}
    \item{\code{GA_EASTWARD_SEPARATION}}{The distance between the center of the sample volume (or the inlet of the intake tube) of the current gas analyzer and the sonic anemometer, as measured horizontallu along the east-west axis (see EddyPro manual).)}
    \item{\code{GA_VERTICAL_SEPARATION}}{The distance between the center of the sample volume (or the inlet of the intake tube) of the current gas analyzer and the sonic anemometer, as measured vertically (see EddyPro manual).}
    \item{\code{GA_TUBE_DIAMETER}}{The inside diameter of the intake tube in centimeters. Mandatory for closed path system.}
    \item{\code{GA_FLOWRATE}}{The flow rate in the intake tube. Mandatory for closed path system.}
    \item{\code{GA_TUBE_LENGTH}}{The length of the intake tube in centimeters.}
    \item{\code{FILE_DURATION}}{File lenght duration in minutes}
    \item{\code{ACQUISITION_FREQUENCY}}{Number of sample records per second.}
    \item{\code{FILE_FORMAT}}{ASCII (fixed)}
    \item{\code{FILE_EXTENSION}}{csv (fixed)}
    \item{\code{LN}}{99 (fixed)}
    \item{\code{FN}}{1 (fixed)}
    \item{\code{EXTERNAL_TIMESTAMP}}{END (fixed)}
    \item{\code{EOL}}{crlf (fixed)}
    \item{\code{SEPARATOR}}{comma (fixed)}
    \item{\code{MISSING_DATA_STRING}}{-9999 (fixed)}
    \item{\code{NROW_HEADER}}{1 (fixed)}
    \item{\code{UVW_UNITS}}{Units of wind speed components (see EddyPro manual).}
    \item{\code{T_SONIC_UNITS}}{Units of sonic temperature (see EddyPro manual).}
    \item{\code{T_CELL_UNITS}}{Units of closed path gas analyzer cell temperature (see EddyPro manual).}
    \item{\code{P_CELL_UNITS}}{Units of closed path gas analyzer cell pressure (see EddyPro manual).}
    \item{\code{CO2_measure_type}}{Carbon dioxide concentration measurement type (e.g. mixing_ratio, see EddyPro manual).}
    \item{\code{CO2_UNITS}}{Units of carbon dioxide concentration (e.g. ppm, see EddyPro manual).}
    \item{\code{H2O_measure_type}}{Water vapor concentration measurement type (e.g. mixing_ratio, see EddyPro manual).}
    \item{\code{H2O_UNITS}}{Units of water vapor concentration (e.g. ppt, see EddyPro manual).}
    \item{\code{SA_DIAG}}{Units of the sonic anemometer diagnostics (dimensionless, fixed)}
    \item{\code{GA_DIAG}}{Units of the gas analyzer diagnostics (dimensionless (fixed).}
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
Metadata elaborated by ICOS-ETC team (Integrated Carbon Observation System (ICOS-ETC) European Research Infrastructure - Ecosystem Thematic Center, \url{http://www.icos-etc.eu/icos/}).
}
\references{
LI-COR Biosciences: EddyPro 7.0.4: Help and User's Guide, LI-COR Biosciences, Lincoln, Nebraska USA, \url{www.licor.com/EddyPro}, 2019.
}
\examples{
data(ecmd_table)
str(ecmd_table)
}
\keyword{datasets}
