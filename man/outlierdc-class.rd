\name{OutlierDC-class}
\Rdversion{1.1}
\docType{class}
\alias{OutlierDC-class}
\title{ "OutlierDC" class}
\description{
	The \code{S4} class of \code{OutlierDC}
}
\section{Objects from the "OutlierDC" Class}{
    Objects can be created by calls of the form \code{new("OutlierDC")}.
}
\section{Slots}{
  \describe{
    \item{\code{call}:}{evaluated function call}
    \item{\code{formula}:}{formula to be used with the type of \code{"Formula"} }
    \item{\code{raw.data}:}{data to be used with the type of \code{"data.frame"} }
    \item{\code{refined.data}:}{data set after removing outliers}
    \item{\code{outlier.data}:}{data set containing outliers}
    \item{\code{coefficients}:}{estimated censored quantile regression coefficient matrix}
    \item{\code{fitted.mat}:}{censored quantile regression fitted value matrix with the type of \code{"matrix"} }
    \item{\code{score}:}{outlying scores (scoring algorithm) or residuals (residual-based algorithm)}
    \item{\code{score.boot}:}{bootstrapping estimation for the outlying scores}
    \item{\code{boot.index}:}{estiamted bootstrap samples used for the outlying scores}
    \item{\code{boot.dist}:}{Empirical quantiles by the Jackknife-after-Bootstrap (JaB) estimation for outlying scores}
    \item{\code{cutoff}:}{estimated scale parameter for the residual-based algorithm}
    \item{\code{lower}:}{lower fence vector used for the boxplot and scoring algorithms with the type of \code{"vector"} }
    \item{\code{upper}:}{upper fence vector used for the boxplot and scoring algorithms with the type of \code{"vector"} }
    \item{\code{outliers}:}{logical vector to determine which observations are outliers}
    \item{\code{n.outliers}:}{number of outliers to be used. The object of class \code{"integer"}. }
    \item{\code{alg}:}{an outlier detection algorithm to be used}
    \item{\code{reg}:}{regression method such as Cox PH or censored quantile regression to be used}
    \item{\code{kr}:}{value to be used for the tightness of cut-offs in the residual-based algorithm}
    \item{\code{kb}:}{numeric value to be used for the tightness of cut-offs in the boxplot algorithm}
    \item{\code{ks}:}{numeric value to be used for the tightness of upper fence cut-offs used for the scoring algorithm with \code{update} function}
    \item{\code{fence}:}{type of fence to be used in the model fitting}
    \item{\code{alpha}:}{numeric value of the significance}
    %\item{\code{LB}:}{lower fence used for the scoring algorithm with \code{update} function. The object of class \code{"numeric"} }
    }
}
\section{Methods}{
  \describe{
    \item{coef}{\code{signature(object = "OutlierDC")}: Print the coefficient matrix of censored quantile regression to be used. See \code{\link{coef}}. }
    \item{plot}{\code{signature(x = "OutlierDC", y = "missing")}: See \code{\link{plot}}. }
    \item{show}{\code{signature(object = "OutlierDC")}: See \code{\link{show}}. }
    \item{update}{\code{signature(object = "OutlierDC")}: Update the fitted object to find outliers in scoring algorithm. See \code{\link{update}}. }
    %\item{jab}{\code{signature(object = "OutlierDC")}: Find outliers using Jackknife-after-Bootstrap sampling for the scoring algorithm. See \code{\link{jab}}. }
	 }
}
\seealso{
    \code{\link{OutlierDC-package}} \cr 
    \code{\link{odc}}, \code{\link{coef}}, \code{\link{plot}}, \code{\link{show}}, \code{\link{update}}, \code{\link{JaB}}
}
\examples{
    showClass("OutlierDC")
}
\keyword{class}
%end