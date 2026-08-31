\name{extract_metrics}
\alias{extract_metrics}
%- Crop point clouds with polygons
\title{
Exrtact features from different Remote Sensing layers.
}
\description{
Extract statistical features from  different kinds of layers, ranging from Raster layers to Point Clouds.
The function is designed to work withone one hand, orthomosaics generated from a DJI P4M drone, but also works with another 5 band raster, where the channels one to five are blue, green, red, Red Edge and Near Infrared.
Here, also texture metrics can be computed, if an 8-band texture layer (simple Haralick features) is passed on. On the other hand, features can be extracted from a point cloud.
Note that here a maximum height should be set to reduce the number of points for extraction.
Paralellisation is can not be turned off for texture and ortho extraction. The extraction of pixelvalues runs with the terra-package and cannot be run parallel.
However, the computation of metrics for each segment based on the pixelwise extract from terra, is paralellised. This increases the speed especially for large umber (>500) of segments.
For smaller number of segments, the paralellisation is not speeding up the process.
}
\usage{
extract_metrics(data_path, type=NULL, date = NULL, polygons=NULL, stepsize = 0.1, out_path = NULL, PC_type = NULL, n_cores = NULL, parallel = TRUE)
}

\arguments{
  \item{data_path}{

} \item{type}{
AAttention: when using LAS, the function assumes that each segment to extract features from is a clipped, seperate LAS file
} \item{date}{

}\item{polygons}{

}\item{stepsize}{

}
}\item{out_path}{

}
}\item{PC_type}{
When extracting from a point cloud (PC), additional specification of the PC format is needed. Use PC_type = 1 for '.txt' format, PC_type = 2 for a colored LiDAR-PC and PC_type = 3 for a multispectral photogrammetric PC (containing a 4th color)
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
