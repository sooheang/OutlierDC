\name{odc}
\alias{odc}
\title{Outlier detection for the analysis of censored data}
\description{
	Detecting outlying observations for censored data, especially with an application of lifetime studies.
}
\usage{
 odc(formula, data, 
  alg = c("score", "boxplot","residual"), 
  reg = c("Wang", "PengHuang", "Portnoy"), 
  fence = c("UB", "LB", "both"),
  kr = 3, kb = 1.5, h = .05)
}
\arguments{
  \item{formula}{
    type of \code{Formula} object with a \code{survival} object on the left-hand side of the ~ operator and covariate terms on the right-hand side. The survival object with survival time and its censoring status is constructed by the \code{\link[survival]{Surv}} function in \code{survival} package.
  }
  \item{data}{
    data frame with variables used in the \code{formula}. It needs at least three variables, including survival time, censoring status, and covariates.
  }
  \item{alg}{
    type of an outlier detection algorithm. Three algorithms are provided. The options \code{"score"}, \code{"boxplot"}, and \code{"residual"} implement the scoring, boxplot, and residual-based algorithms, respectively. The default algorithm is the scoring algorithm with an option \code{"score"}.

  }
  \item{reg}{
    type of a regression method used as a basis for outlier detection algorithms. The options \code{"Wang"}, \code{"Portnoy"}, and \code{"PengHuang"} conduct Wang and Wang's, Portnoy's, and Peng and Huang's censored quantile regression approaches, respectively. The default is a local CQR with an option \code{"Wang"}.
  }
  \item{fence}{
    type of an outlying fence. Three options are provided. The option \code{"UB"} provides the upper fence to find an observation who "lived far too long". The option \code{"LB"} provides the lower fence to find an observation who "died far too early". The option \code{"both"} provides the upper and lower fences, simultaneously. 
  }
  \item{kr}{
    numeric value to control the tightness of cut-offs for the residual algorithm with a default value of 3.
  }
  \item{kb}{
    numeric value to control the tightness of cut-offs for the boxplot algorithm with a default value of 1.5.
  }
   \item{h}{
   bandwidth for locally weighted censored quantile regression with a default value of 0.05.
  }
}
\details{
	The \code{odc} function conducts three outlier detection algorithms on the basis of censored quantile regression. Three outlier detection algorithms were implemented: residual-based, boxplot, and scoring algorithms. The residual-based algorithm detects outlying observations using constant scale estimates; however, it does not account for the heterogeneity of variability. When the data is extremely heterogeneous, the boxplot algorithm with censored quantile regression is more effective. The residual-based and boxplot algorithms produce cut-offs to determine whether observations are outliers. In contrast, the scoring algorithm provides the outlying magnitude or deviation of each point from the distribution of observations. Outlier detection is achieved by visualising the scores. 
}
\value{
	an object of the S4 class "OutlierDC" with the following slots:\cr
	call: evaluated function call\cr
  formula: formula to be used\cr
  raw.data: data to be used for model fitting\cr
  refined.data: the data set after removing outliers\cr
  refined.data: the data set containing outliers\cr
  coefficients: the estimated censored quantile regression coefficient matrix consisting of 10th, 25th, 50th, 75th, and 90th quantiles\cr
  fitted.mat: the censored quantile regression fitted value matrix consisting of 10th, 25th, 50th, 75th, and 90th quantiles\cr 
  score: outlying scores (scoring algorithm) or residuals (residual-based algorithm)\cr
  cutoff: estimated scale parameter for the residual-based algorithm\cr
  lower: lower fence vector used for the boxplot and scoring algorithms\cr
  upper: upper fence vector used for the boxplot and scoring algorithms\cr
  outliers: logical vector to determine which observations are outliers\cr
  n.outliers: number of outliers detected\cr
  method: outlier detection method to be used\cr
  rq.model: censored quantile regression to be used\cr
  kr: a value to be used for the tightness of cut-offs in the residual algorithm \cr
  kb: a value to be used for the tightness of cut-offs in the boxplot algorithm \cr
  ks: a value to be used for the tightness of upper fence cut-offs used for the scoring algorithm with the \code{\link{update}} function\cr
  fence: type of fence to be used in the model fitting \cr
  alpha: numeric value for the significance level \cr
  boot.dist: empirical quantiles by Jackknife-after-Bootstrapping 
  %B: a value to be used for the scoring algorithm in order to determine the empirical distribution of each observation.
  %LB: lower fence used for the scoring algorithm with the \code{\link{update}} function\cr
}
\source{
Eo, S-H, Hong, S-M, and Cho, H. (2014+). Identification of outlying observations with quantile regression for censored data, \emph{Submitted}.

Martin, M. A., and Roberts, S. (2010). Jackknife-after-bootstrap regression influence diagnostics, \emph{Journal of Nonparametric Statistics}, 22, 257-269.

Nardi, A., and Schemper, M. (1999). New residuals for Cox regression and their application to outlier screening, \emph{Biometrics}, 55, 523-529.

Wang H. J., and Wang, L. (2009). Locally weighted censored quantile regression. \emph{Journal of the American Statistical Association}, 104, 1117-1128. 
}
\seealso{
    \code{\link{OutlierDC-package}} \cr 
    \code{\link{coef}}, \code{\link{plot}}, \code{\link{show}}, \code{\link{update}}
}
\examples{
  \dontrun{

library(OutlierDC)
data(ebd)
str(ebd)

####
# outlier detection using the scoring algorithm
fit = odc(Surv(log(time), status) ~ meta, data = ebd)
fit

# A threshold is added by k_s to this plot using the updata() function
fit1 = update(fit, ks = 4)
fit1
plot(fit1)

# A threshold can be determined by using the empirical distribution of standard deviation for the outlying scores
fit2 = bodc(fit, B = 500)
fit2

####
# outlier detection using the residual-based algorithm
fit3 = odc(Surv(log(time), status) ~ meta, data = ebd, alg = "residual")
fit3
plot(fit3, main = "Residual-based algorithm")


fit3@outlier.data

####
# outlier detection using the boxplot-based algorithm
fit4 = odc(Surv(log(time), status) ~ meta, data = ebd, alg = "boxplot")
fit4
plot(fit4, main = "Boxplot-based algorithm", xlab = "Number of metastatic lymph nodes",, ylab = "Log of survival times")
}
}
\keyword{odc}


