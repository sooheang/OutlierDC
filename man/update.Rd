\name{update}
\docType{methods}
\alias{update}
\alias{update,OutlierDC-method}

\title{ Update cut-off values for scoring algorithm by controlling the absolute value of upper and lower fences or the significane level alpha}

\description{
This function updates a cut-off values for scoring algorithm by controlling the absolute value of upper and lower fences, or the significane level alpha for Jackknife-after-Bootstrapping (JaB) distribution. Using the call stored in the object, the \code{update} function declares outlying observatoins based on the QQ plot. \code{alpha} is used to control the significant level of the JaB distribution, and \code{ks} and \code{LB} are used to set the upper and lower fences, respectively.
% and \code{LB} is used to set the lower cut-off bound. 
}
\usage{
    \S4method{update}{OutlierDC}(object, alpha = 0.01, ks = NA, LB = NA)
}
\arguments{
    \item{object}{
        fitted model object of class \code{\linkS4class{OutlierDC}}.
    }
    \item{alpha}{
        significant level for the Jackknife-after-Bootstrap resamplings. This arguments is only VALID after running a \code{\link{JaB}()} function. 
    }
    \item{ks}{ 
		cut-off value for the upper fence.
	}
    \item{LB}{  
		cut-off value for the lower fence.
	}
}
\details{
  This function is a generic function called \code{update} for the S4 class \code{OutlierDC}. Cut-off bounds are added to find outliers on the QQ plot. 
}
\seealso{
    \code{\link{odc}} and \code{\linkS4class{OutlierDC}} class
}
\keyword{methods}
