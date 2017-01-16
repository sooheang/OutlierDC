\name{coef}
\docType{methods}
\alias{coef}
\alias{coef,OutlierDC-method}
\title{ a \code{coef} method for "OutlierDC".}
\description{
    \code{coef} is a generic function which extracts model coefficient matrix including the 5th, 10th, 25th, 50th, 75th, 90th, and 95th quantile estimates. 
}
\usage{
    \S4method{coef}{OutlierDC}(object)
}
\arguments{
    \item{object}{
        an object with class \code{\linkS4class{OutlierDC}}.  
    }
}
\details{
  This function is a generic function \code{coef} for the S4 class \code{OutlierDC}.  It can be invoked by calling print for an object of the appropriate class, or directly by calling \code{coef} regardless of the class of the object. 
}
\seealso{
    \code{\link{odc}} and \code{\linkS4class{OutlierDC}} class
}
\keyword{methods} 