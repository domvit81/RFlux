\name{RFlux-package}
\alias{RFlux-package}
\alias{RFlux}
\docType{package}
\title{
\packageTitle{RFlux}
}
\description{
\packageDescription{RFlux}
}
\details{

RFlux package ingests eddy covariance rawdata sampled by either open- or closed-path system and implement the processing pipeline adopted by the ICOS-ETC (Integrated Carbon Observation System European Research Infrastructure - Ecosystem Thematic Center, \url{http://www.icos-etc.eu/icos/}). All metadata information have to be reported in the \emph{filename}_ecmd.csv table. Such information are then processed by the \code{\link{get_md}} function.

The processing of rawdata aims at
\describe{
\item{i}{estimating fluxes and other micrometeorolgical parameters.}
\item{ii}{performing data quality control.}
}

Flux estimation involves the following options/methods:
\itemize{
\item Planar fit for tilt correction (Wilczak et al, 2001).
\item Maximum cross-covariance method for time lag determination (see Rebmann et al, 2012).
\item Block average (see Moncrieff et al, 2004).
\item WPL correction, only for open path systems (Webb et al, 1980).
\item In-situ spectral corrections (Fratini et al, 2012).
}

The open source EddyPro software (LI-COR Biosciences, 2019) is used to this aim employing also the estimation of micrometeorological parameters useful in subsequent analyses. It is required the EddyPro software is installed on your system (for download see \url{www.licor.com/EddyPro}).

Quality control involves the data cleaning procedure described in Vitale et al (2019). Its implementation involves a three-step procedure
\describe{
\item{Step 1:}{Estimation of the test statistics via the \code{\link{qcStat}} function.}
\item{Step 2:}{Generating the workset via the \code{\link{ecworkset}} function.}
\item{Step 3:}{Application of data cleaning procedure (including despiking) via the \code{\link{cleanFlux}} function.}
}

\emph{Acknowledgements.} Rflux has been developed in the context of the ICOS Ecosystem Thematic Centre. DV thanks the ENVRIPLUS H2020 European project (Grant Agreement 654182) for the support. DP thanks the ENVRIFAIR H2020 European project (Grant Agreement 824068) for the support.



%%The DESCRIPTION file:
%%\packageDESCRIPTION{RFlux}
%%\packageIndices{RFlux}
%%~~ An overview of how to use the package, including the most important ~~
%%~~ functions ~~
}
\author{
\packageAuthor{RFlux}

Maintainer: \packageMaintainer{RFlux}
}
\references{
Fratini, G., Ibrom, A., Arriga, N., Burba, G., Papale, D. (2012). Relative humidity effects on water vapour fluxes measured with closed-path eddy-covariance systems with short sampling lines. Agricultural and forest meteorology, 165, pp 53-63, doi: \url{https://doi.org/10.1016/j.agrformet.2012.05.018}.

LI-COR Biosciences: EddyPro 7.0.4: Help and User's Guide, LI-COR Biosciences, Lincoln, Nebraska USA, \url{www.licor.com/EddyPro}, 2019.

Moncrieff, J., Clement, R., Finnigan, J., Meyers, T. (2004). Averaging, detrending, and filtering of eddy covariance time series. In Handbook of micrometeorology, pp. 7-31, Springer, Dordrecht, doi: \url{https://doi.org/10.1007/1-4020-2265-4_2}.

Rebmann, C., Kolle, O., Heinesch, B., Queck, R., Ibrom, A., Aubinet, M. (2012). Data acquisition and flux calculations. In Eddy covariance, pp. 59-83, Springer, Dordrecht, doi: \url{https://doi.org/10.1007/978-94-007-2351-1_3}.

Vitale, D. Fratini, G. Bilancia, M. Nicolini, G. Sabbatini, S. Papale, D. (2019).
A robust data cleaning procedure for eddy covariance flux measurements,
Biogeosciences Discussions, pp 1-36, doi: \url{https://doi.org/10.5194/bg-2019-270}.

Webb, E.K., Pearman, G.I., Leuning, R. (1980). Correction of flux measurements for density effects due to heat and water vapour transfer. Quarterly Journal of the Royal Meteorological Society, 106(447), pp 85-100, doi: \url{https://doi.org/10.1002/qj.49710644707}.

Wilczak, J.M., Oncley, S.P., Stage, S.A. (2001). Sonic anemometer tilt correction algorithms. Boundary-Layer Meteorology, 99(1), pp 127-150, doi: \url{https://doi.org/10.1023/A:1018966204465}.
}
%%~~ Optionally other standard keywords, one per line, from file KEYWORDS ~~
%%~~ in the R documentation directory ~~
\keyword{ package }
\seealso{
%%~~ Optional links to other man pages, e.g. ~~
%%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
}
\examples{
%%~~ simple examples of the most important functions ~~
}
