\name{rad_Cor_P4M}
\alias{Radiometric_Correction_P4M}
%- Crop point clouds with polygons
\title{
Correct Orthomosaics from a P4M UAV with given correction factors from a spctroradiometer
}
\description{
This function extracts the average values of polygons that represent reflactance calibration targets in the field. Coefficients based on a linear model are extracted and applied as correction factors to the input raster.
}
\usage{
rad_Cor_P4M(spath = NULL, rpath = NULL, spectral_file = NULL)
}

\arguments{
  \item{spath}{
Path to shape file covering the inner pixels of a reference target on the uncorrected orthomosaic. Needs to contain a field "Target_ID" with a unique ID and "Farbe" describing the color of the target.
} \item{rpath}{
Path to processed multiband orthomosaic containing the five spectral bands of the DJI P4M.
} \item{spectral_file}{
Path of the .csv file containing the reference measures from the spctroradiometer.

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
The function will write a radiometrically corrected tif-file to the input directory of the uncorrected raster.
}
\references{
%% ~put references to the literature/web site here ~
}\author{

Author: MrFlywheel
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
