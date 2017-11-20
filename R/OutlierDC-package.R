#' Functions for detecting outlying observations in the analysis of censored
#' data, especially with an application to lifetime studies.
#' 
#' This package offers exploratory outlier detection algorithms for the
#' analysis of censored data, especially, with an application to lifetime
#' studies.
#' 
#' \tabular{ll}{ Package: \tab OutlierDC\cr Type: \tab Package\cr Version: \tab
#' 1.0.0\cr Date: \tab 2015-06-06\cr License: \tab GPL version 3\cr LazyLoad:
#' \tab no\cr }
#' 
#' @name OutlierDC-package
#' @aliases OutlierDC OutlierDC-package
#' @docType package
#' @note We would like to thank Huxia Judy Wang and Lan Wang for the permission
#' to use the function LCRQ().
#' @author Soo-Heang EO and HyungJun CHO \cr Maintainer: Soo-Heang EO
#' <eo.sooheang@@gmail.com>
#' @seealso \code{\link{odc}}, \code{\link{plot}}, \code{\link{coef}},
#' \code{\link{show}}, \code{\link[quantreg:crq]{quantreg}}
#' @references Eo, S-H, Hong, S-M, and Cho, H. (2015+). Identification of
#' outlying observations with quantile regression for censored data,
#' \emph{Submitted}.
#' 
#' Balakrishnan, N., Chimitova, E., Galanova, N., & Vedernikova, M. (2013).
#' Testing Goodness of Fit of Parametric AFT and PH Models with Residuals.
#' \emph{Communications in Statistics - Simulation and Computation}, 42(6),
#' 1352–1367.
#' 
#' Wang H. J., and Wang, L. (2009). Locally weighted censored quantile
#' regression. \emph{Journal of the American Statistical Association}, 104,
#' 1117-1128.
#' 
#' Martin, M. A., and Roberts, S. (2010). Jackknife-after-bootstrap regression
#' influence diagnostics, \emph{Journal of Nonparametric Statistics}, 22,
#' 257-269.
#' 
#' Nardi, A., and Schemper, M. (1999). New residuals for cox regression and
#' their application to outlier screening, \emph{Biometrics}, 55, 523-529.
#' @keywords package
NULL





#' a \code{coef} method for "OutlierDC".
#' 
#' \code{coef} is a generic function which extracts model coefficient matrix
#' including the 5th, 10th, 25th, 50th, 75th, 90th, and 95th quantile
#' estimates.
#' 
#' This function is a generic function \code{coef} for the S4 class
#' \code{OutlierDC}.  It can be invoked by calling print for an object of the
#' appropriate class, or directly by calling \code{coef} regardless of the
#' class of the object.
#' 
#' @name coef
#' @aliases coef coef,OutlierDC-method
#' @docType methods
#' @param object an object with class \code{\linkS4class{OutlierDC}}.
#' @seealso \code{\link{odc}} and \code{\linkS4class{OutlierDC}} class
#' @keywords methods
NULL





#' Extrahepatic Cholangiocarcinoma Data
#' 
#' The extrahepatic cholangiocarcinoma data comes form the US Surveilance,
#' Epidemiology, and End Results (SEER) program of the National Cancer
#' Institute.
#' 
#' 
#' @name ebd
#' @docType data
#' @source Hankey B., Ries L., and Edwards B. (1999). The surveillance,
#' epidemiology, and end results program a national resource. \emph{Cencer
#' Epidemiology Biomarkers and Prevention}, 12:1117-1121.
#' @keywords datasets
#' @examples
#' 
#' 	data(ebd)
#' 
NULL





#' "OutlierDC" class
#' 
#' The \code{S4} class of \code{OutlierDC}
#' 
#' 
#' @name OutlierDC-class
#' @docType class
#' @section Objects from the "OutlierDC" Class: Objects can be created by calls
#' of the form \code{new("OutlierDC")}.
#' @seealso \code{\link{OutlierDC-package}} \cr \code{\link{odc}},
#' \code{\link{coef}}, \code{\link{plot}}, \code{\link{show}},
#' \code{\link{update}}, \code{\link{JaB}}
#' @keywords class
#' @examples
#' 
#'     showClass("OutlierDC")
#' 
NULL





#' a plot-method for a "OutlierDC" object
#' 
#' This function provides three different results. If the algorithm is
#' \code{score}, it draws a normal quantile-quantile plot of the outlying
#' scores. If the algorithm is \code{boxplot}, the scatter plot of log survival
#' times against the covariate used is given. Lastly, if the algorithm is
#' \code{residual}, it offers a residual plot.
#' 
#' 
#' @name plot
#' @aliases plot plot,OutlierDC-method
#' @docType methods
#' @param x fitted model object of class \code{\linkS4class{OutlierDC}}.
#' @param y missing value. Not used parameter in this package.
#' @param ...  \code{\link[graphics:plot.default]{plot.default}} arguments
#' @seealso \code{\link{odc}} and \code{\linkS4class{OutlierDC}} class
#' @keywords methods
NULL





#' a show method for \code{OutlierDC}
#' 
#' This function provides a summary for the \code{OutlierDC} class.
#' 
#' This function is a method for the generic function \code{show} for the S4
#' class \code{OutlierDC}.  It can be invoked by calling print for an object of
#' the appropriate class, or directly by calling \code{show} regardless of the
#' class of the object.
#' 
#' @name show
#' @aliases show show,OutlierDC-method
#' @docType methods
#' @param object fitted object of class \code{\linkS4class{OutlierDC}}.
#' @seealso \code{\link{odc}} and \code{\linkS4class{OutlierDC}}
#' @keywords methods
NULL





#' Toy dataset
#' 
#' a \code{toy} dataset
#' 
#' 
#' @name toy
#' @docType data
#' @return \item{time}{time to event} \item{status}{censoring status}
#' \item{meta}{a covariate variable that describes the number of metastatic
#' lymph nodes}
#' @keywords datasets
NULL





#' Update cut-off values for scoring algorithm by controlling the absolute
#' value of upper and lower fences or the significane level alpha
#' 
#' This function updates a cut-off values for scoring algorithm by controlling
#' the absolute value of upper and lower fences, or the significane level alpha
#' for Jackknife-after-Bootstrapping (JaB) distribution. Using the call stored
#' in the object, the \code{update} function declares outlying observatoins
#' based on the QQ plot. \code{alpha} is used to control the significant level
#' of the JaB distribution, and \code{ks} and \code{LB} are used to set the
#' upper and lower fences, respectively. % and \code{LB} is used to set the
#' lower cut-off bound.
#' 
#' This function is a generic function called \code{update} for the S4 class
#' \code{OutlierDC}. Cut-off bounds are added to find outliers on the QQ plot.
#' 
#' @name update
#' @aliases update update,OutlierDC-method
#' @docType methods
#' @param object fitted model object of class \code{\linkS4class{OutlierDC}}.
#' @param alpha significant level for the Jackknife-after-Bootstrap
#' resamplings. This arguments is only VALID after running a
#' \code{\link{JaB}()} function.
#' @param ks cut-off value for the upper fence.
#' @param LB cut-off value for the lower fence.
#' @seealso \code{\link{odc}} and \code{\linkS4class{OutlierDC}} class
#' @keywords methods
NULL



