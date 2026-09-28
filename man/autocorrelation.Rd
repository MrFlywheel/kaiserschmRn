\name{autocorrelation}
\alias{autocorrelation}
%- Crop point clouds with polygons
\title{
Get the estimated autocorrelation range for a given raster.
}
\description{
As many vectors containing integer values can be passed on to search the intersection of all vectors, which will be returned.
}
\usage{
autocorrelation(rast, sample_size = 10000, cut_off = 50, width = 1)
}

\arguments{
  \item{rast}{
One raster layer where autocorrelation range should be estimated
}\item{sample_size}{
The number of pixel values to be sampled. Default is 10.000
}\item{cut_off}{
Cut_off parameter (i.e. range of the x-Axis) to be ploted.
}\item{width}{
Width of bins to compare semivariance estimates
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
The function returns a plot of the semi-variogramm, where the range can be inferred.
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
