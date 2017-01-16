\name{JaB}
\alias{JaB}
\title{Jackknife-after-Bootstrap resampling for scoring algorithm}
\description{
	This function constructs cut-off values for the scoring algorithm using Jankknife-after-Bootstrap resampling.
}
\usage{
 JaB(object, B = 1000, alpha = 0.05, fence = c("UB", "LB", "both"))
}
\arguments{
  \item{object}{
    fitted object of class \code{\linkS4class{OutlierDC}}  
  }
  \item{B}{
    number of Jackknife-after-Bootstrap resamples for the scoring algorithm. The default value is 300
  }
  \item{alpha}{
    magnitude of significance
  }
  \item{fence}{
    a character value to control an outlying fence. three options are provided. The option \code{"UB"} provides the upper fence to find an observation who "lived far too long". The option \code{"LB"} provides the lower fence to find an observation who "died far too early". The option \code{"both"} provides the upper and lower fences simultaneously. 
  }
}
\details{
	The function \code{JaB} provides cut-off values for scoring algorithm using Jackknife-after-Bootstrap resampling techniques.
}
\seealso{
    \code{\link{OutlierDC-package}} \cr 
    \code{\link{odc}}, \code{\link{bodc}} \cr
    \code{\link{coef}}, \code{\link{plot}}, \code{\link{show}}, \code{\link{summary}}, \code{\link{update}}
}
\keyword{odc}


