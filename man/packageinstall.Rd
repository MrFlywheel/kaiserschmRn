\name{package_install}
\alias{package_install}
\title{Check if packages are installed. If not, do so.}
\usage{
package.install(x)
}
\arguments{
\item{x}{Package}
}
\description{
Check if packages are installed. If not, do so.
}
\examples{
library(kaiserschmRn)

# check/install one package
package.install("kaiserschmRn")

# check/install multiple packages
package.install(c("kaiserschmRn","kaiserschmRn"))
}
