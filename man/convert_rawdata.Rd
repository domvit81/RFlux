\name{convert_rawdata}
\alias{convert_rawdata}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Raw Data File Conversion}
\description{
Converting a generic ASCII raw data file in a format suitable for RFlux package
}
\usage{
convert_rawdata(
file_path_in,
file_path_out, 
info_U=list(NULL,NULL,1,0), 
info_V=list(NULL,NULL,1,0),
info_W=list(NULL,NULL,1,0),
info_T_SONIC=list(NULL,NULL,1,0),
info_CO2=list(NULL,NULL,1,0),
info_H2O=list(NULL,NULL,1,0),
info_T_CELL=list(NULL,NULL,1,0),
info_T_CELL_IN=list(NULL,NULL,1,0),
info_T_CELL_OUT=list(NULL,NULL,1,0), 
info_PRESS_CELL=list(NULL,NULL,1,0), 
info_SA_DIAG=NULL, 
info_GA_DIAG=NULL, 
na.strings, 
nrow.header, 
sep, 
dec, 
timestamp_loc=c(NULL, NULL), 
timestamp_format, 
siteID="Cc-Xxx", 
labels=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file_path_in}{
  path where raw data file is stored.
  }
  \item{file_path_out}{
  path where converted raw data file will be stored.
	}
	\item{info_U}{
	a list of 4 elements denoting the column position (integer) of the horizontal U wind component, its units (character string, see details), the multiplier and the offset (1 and 0, respectively, in case no linear conversion should be applied).
	}
	\item{info_V}{
	a list of 4 elements denoting the column position (integer) of the horizontal V wind component, its units (character string, see details), the multiplier and the offset (1 and 0, respectively, in case no linear conversion should be applied).
	}
	\item{info_W}{
	a list of 4 elements denoting the column position (integer) of the vertical W wind component, its units (character string, see details), the multiplier and the offset (1 and 0, respectively, in case no linear conversion should be applied).
	}
	\item{info_T_SONIC}{
	a list of 4 elements denoting the column position (integer) of the sonic temperature, its units (character string, see details), the multiplier and the offset (1 and 0, respectively, in case no linear conversion should be applied).
	}
	\item{info_CO2}{
	a list of 4 elements denoting the column position (integer) of CO2 concentration, its units (character string, see details), the multiplier and the offset (1 and 0, respectively, in case no linear conversion should be applied).
	}
	\item{info_H2O}{
	a list of 4 elements denoting the column position (integer) of the H2O concentration, its units (character string, see details), the multiplier and the offset (1 and 0, respectively, in case no linear conversion should be applied).
	}
	\item{info_T_CELL}{
	a list of 4 elements denoting the column position (integer) of cell temperature (only for EC closed path systems), its units (character string, see details), the multiplier and the offset (1 and 0, respectively, in case no linear conversion should be applied).
	}
	\item{info_T_CELL_IN}{
	a list of 4 elements denoting the column position (integer) of in-cell temperature (only for EC closed path systems), its units (character string, see details), the multiplier and the offset (1 and 0, respectively, in case no linear conversion should be applied).
	}
	\item{info_T_CELL_OUT}{
	a list of 4 elements denoting the column position (integer) of out-cell temperature (only for EC closed path systems), its units (character string, see details), the multiplier and the offset (1 and 0, respectively, in case no linear conversion should be applied).
	}
	\item{info_PRESS_CELL}{
	a list of 4 elements denoting the column position (integer) of cell pressure (only for EC closed path systems), its units (character string, see details), the multiplier and the offset (1 and 0, respectively, in case no linear conversion should be applied).
	}
	\item{info_SA_DIAG}{
	integer. Column position of the sonic anemometer diagnostics.
	}
	\item{info_GA_DIAG}{
	integer. Column position of the gas analyzer diagnostics.
	}
	\item{na.strings}{
	a character string denoting missing values in the raw data file.
	}
	\item{nrow.header}{
	the number of header rows in the raw data file.
	}
	\item{sep}{
	the field separator character of the raw data file. See \code{\link{read.table}} for specifications.
	}
	\item{dec}{
	the character used in the raw data file for decimal points.
	}
	\item{timestamp_loc}{
	a vector of 2 negative integers denoting the location of the timestamp in the raw data file name starting from the right (eg. -21, -2).
	}
	\item{timestamp_format}{
	a character string denoting the format of the timestamp included in the raw data file name. It must be specified according to specifications of the \code{\link{strptime}} function.)
	}
	\item{siteID}{
	a 6-character string denoting the Site ID (eg. CC-Xxx).
	}
	\item{labels}{
	(optional) a character string denoting further information. Default is NULL.
	}
}
\details{
Possible units specification in the info_XX field are one among:
mm s-1, cm s-1, m s-1 for U, V, W wind components;
K, C for sonic temperature (T_SONIC) and cell temperatures (T_CELL, T_CELL_IN, T_CELL_OUT);
ppt, ppm, ppb, umol m-3, mmol m-3, ug m-3, mg m-3, g m-3 for CO2, H2O concentration variables;
Pa, hPa, KPa for cell pressure (PRESS_CELL).
The multiplier and the offset values in info_XX field aim at converting data expressed in Volt in their physical units. They must be set equals to 1 and 0, respectively, in case no (linear) conversion is required.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
A .csv file stored in the file_path_out containing the following variables: U (m s-1), V (m s-1), W (m s-1), T_SONIC (K), CO2 (mmol m-3 in case of molar density, umol/mol (ppm) in case of molar fraction/mixing ratio), H2O (mmol m-3 in case of molar density, mmol/mol (ppt) in case of molar fraction/mixing ratio), SA_DIAG (dimensionless), GA_DIAG (dimensionless), T_CELL (C), T_CELL_IN (C), T_CELL_OUT (C), PRESS_CELL (KPa). Missing values are denoted by "-9999".
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
##Single raw data file conversion (TO DO)
##Block of raw data file conversion in parallel mode (TO DO)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{datasets}% use one of  RShowDoc("KEYWORDS")
