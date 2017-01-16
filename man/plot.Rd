\name{plot}
\docType{methods}
\alias{plot}
\alias{plot,OutlierDC-method}
\title{a plot-method for a "OutlierDC" object}
\description{
    This function provides three different results. If the algorithm is \code{score}, it draws a normal quantile-quantile plot of the outlying scores. If the algorithm is \code{boxplot}, the scatter plot of log survival times against the covariate used is given. Lastly, if the algorithm is \code{residual}, it offers a residual plot. 
}
\usage{
    \S4method{plot}{OutlierDC}(x, y = NA, ...)
}
\arguments{
    \item{x}{
        fitted model object of class \code{\linkS4class{OutlierDC}}. 
    }
    \item{y}{
      missing value. Not used parameter in this package.
    }
    \item{...}{ \code{\link[graphics:plot.default]{plot.default}} arguments}
}
\seealso{
    \code{\link{odc}} and \code{\linkS4class{OutlierDC}} class
}
\keyword{methods}