\name{inst_prob_test}
\alias{inst_prob_test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Instrumental problem detection 
}
\description{
A set of tests aims at detecting instrumental malfunctions affecting eddy covariance systems.
}
\usage{
inst_prob_test(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  raw high frequency eddy covariance time series.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{Skew}{Skewness.}
  \item{Kurt}{Kurtosis index on original data.}
  \item{KID0}{Kurtosis index on difference data.}
  \item{KID1}{Kurtosis index on difference data excluding low resolution problems.}
  \item{HFx}{Homogeneity test statistic based on fluctuation.}
  \item{HDx}{Homogeneity test statistic based on difference data.}
}
\references{
Vitale, D. Fratini, G. Bilancia, M. Nicolini, G. Sabbatini, S. Papale, D.
A robust data cleaning procedure for eddy covariance flux measurements,
Biogeosciences Discussions, 2019, pp 1-36, doi: \url{https://doi.org/10.5194/bg-2019-270}.
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
data(closed_path_rawdata)
inst_prob_test(closed_path_rawdata$W)
}
\keyword{univar}% use one of  RShowDoc("KEYWORDS")
