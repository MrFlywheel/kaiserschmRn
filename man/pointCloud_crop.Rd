\name{pointCloud_crop}
\alias{pointCloud_crop}
%- Crop point clouds with polygons
\title{
Crop '.las' files with polygons
}
\description{
Function to crop '.las' files to the extent of polygons in a very fast way by parallelisation. Uses the crop_ROI function of the lidR package (CRAN.R-project.org/package=lidR), but runs much faster.
}
\usage{
pointCloud_crop(las_path=NULL, shp_path=NULL, out_path=NULL, n_cores=NULL)
}

\arguments{
  \item{las_path}{
Item to be cropped. A full path to a '.las' file. Will not be read to memory directly.
} \item{shp_path}{
Item used for croppping. A full path to a the shape file used for cropping.
} \item{out_path}{
Path of the output directory. Will create a sub folder 'export' in the 'las_path' directory, if not specified.
}\item{n_cores}{
Number of course used for parallelisation. As a default half of the PC's logical cores are used.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
The function will write '.las' files for each polygon.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
MrFlywheel
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x)
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
