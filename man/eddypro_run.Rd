\name{eddypro_run}
\alias{eddypro_run}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fluxes estimation
}
\description{
Estimates flux values and other micrometeorological parameters through a call to LI-COR EddyPro software.
}
\usage{
eddypro_run(siteID, path_eddypro_bin, path_eddypro_projfiles, showLOG = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{siteID}{
  Character string indicative of the site's ID (i.e. CC-Xxx)  
}
  \item{path_eddypro_bin}{
  path where eddypro_rp and eddypro_fcc executables (i.e. the bin folder) are stored.
}
  \item{path_eddypro_projfiles}{
  path where the EddyPro project files (i.e. the \emph{filename}.eddypro file) is stored.
}
  \item{showLOG}{
  logical. If TRUE (default), the EddyPro Output Console Page is shown on R console.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A set of .csv files (e.g. fulloutput, qcdetails and metadata) generated by EddyPro software will be stored in the path specified by \code{\link{get_md}}. %%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
LI-COR Biosciences: EddyPro 7.0.4: Help and User's Guide, LI-COR Biosciences, Lincoln, Nebraska USA, \url{www.licor.com/EddyPro}, 2019.
}
\author{
Domenico Vitale, Dario Papale}
\note{
}

\section{Warning}{
It is required LI-COR EddyPro software is currently installed on your system (see \url{www.licor.com/EddyPro}).
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
\dontrun{
workdir <- getwd()
siteID <- "DE-HoH"
PATH_ECMD <- system.file("extdata", "DE-HoH_ecmd.csv", package = "RFlux")
PATH_RAWDATA <- system.file("extdata", package = "RFlux")
PATH_OUTPUT <- paste0(workdir, "/eddypro/processing")
PATH_SA_FILE <- system.file("extdata", "spectral_assessment.txt", package = "RFlux")
PATH_PF_FILE <- system.file("extdata", "planar_fit.txt", package = "RFlux")

MD <- get_md(path_ecmd=PATH_ECMD,
             path_rawdata=PATH_RAWDATA,
             path_output=PATH_OUTPUT,
             online=TRUE,
             path_sa_file=PATH_SA_FILE,
             path_pf_file=PATH_PF_FILE,
             tlag_meth=2,
             despike_meth="VM97",
             detrend_meth="BA",
             tilt_correction_meth="PF")

PATH_EDDYPRO_BIN <- "/Applications/eddypro.app/Contents/MacOS/bin" ## put your path!
PATH_EDDYPRO_PROJ <- PATH_OUTPUT
eddypro_run(siteID="DE-HoH",
            path_eddypro_bin=PATH_EDDYPRO_BIN,
            path_eddypro_projfiles=PATH_EDDYPRO_PROJ,
            showLOG=TRUE)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{interface}% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
