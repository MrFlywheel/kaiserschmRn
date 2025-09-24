\name{roughness_RMSE}
\alias{roughness_RMSE}
%- Crop point clouds with polygons
\title{
Compute surface roughness of a point cloud
}
\description{
Surface roughness is computed with the RMSE of of z-vlaues compared to a fitted plane, where 'resolution' defines the resolution of a raster cell where the plane is fitted. Needs to be a normalised point cloud, as otherwise the RMSE would vary depending on the location of the cell. This function normalises a point cloud by using lidR's functions when a dtm is provided.
}
\usage{
roughness_RMSE(las_object, res = 1, dtm = NULL, filter_height = NULL)
}

\arguments{
  \item{las_object}{
Point cloud to be used to cumpute roughness. An objec of type 'las'. To avoid long waiting times in case of some error in the function, this las object should be read into memroy already beforehand.
} \item{resolution}{
Resolution of the raster cells to fit the plane.
} \item{dtm}{
Path to a DTM. Only needs to be provided, if point cloud is not already height normalised.
}\item{filter_height}{
Heigh threshold to consider for roughness calculation.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Will plot and return a raster of the surface roughness
}
\references{
lidR package
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

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
