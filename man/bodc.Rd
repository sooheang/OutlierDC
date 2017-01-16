\name{bodc}
\alias{bodc}
\title{Calculate a distribution-free cut-off value for the scoring algorithm}
\description{
	This function constructs cut-off values for the scoring algorithm using Bootstrap distribution of standard deviation.
}
\usage{
 bodc(object, B = 500, q = 0.05, fence = c("UB", "LB", "both"))
}
\arguments{
  \item{object}{
    fitted object of class \code{\linkS4class{OutlierDC}}  
  }
  \item{B}{
    number of bootstrap resamples for the scoring algorithm. The default value is 500.
  }
  \item{q}{
    numeric value of probability with values in [0,1].
  }
  \item{fence}{
    a character value to control an outlying fence. three options are provided. The option \code{"UB"} provides the upper fence to find an observation who "lived far too long". The option \code{"LB"} provides the lower fence to find an observation who "died far too early". The option \code{"both"} provides the upper and lower fences simultaneously. 
  }
}
\details{
	The function \code{bodc} provides a distribution-free threshold for scoring algorithm by using the empirical distribution of the standard deviation of outlying scores.
}
\seealso{
    \code{\link{OutlierDC-package}} \cr 
    \code{\link{odc}}, \code{\link{JaB}}\cr
    \code{\link{coef}}, \code{\link{plot}}, \code{\link{show}}, \code{\link{summary}}, \code{\link{update}}
}
\keyword{odc}


