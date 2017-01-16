\name{OutlierDC-package}
\Rdversion{1.1}
\docType{package}
\alias{OutlierDC}
\alias{OutlierDC-package}
\title{
	Functions for detecting outlying observations in the analysis of censored data, especially with an application to lifetime studies.
}
\description{
  This package offers exploratory outlier detection algorithms for the analysis of censored data, especially, with an application to lifetime studies.
}
\details{
    \tabular{ll}{
        Package: \tab OutlierDC\cr
        Type: \tab Package\cr
        Version: \tab 1.0.0\cr
        Date: \tab 2015-06-06\cr
        License: \tab GPL version 3\cr
        LazyLoad: \tab no\cr
    }
}
\author{
Soo-Heang EO and HyungJun CHO \cr
Maintainer: Soo-Heang EO <eo.sooheang@gmail.com>
}
\references{
Eo, S-H, Hong, S-M, and Cho, H. (2015+). Identification of outlying observations with quantile regression for censored data, \emph{Submitted}.

Balakrishnan, N., Chimitova, E., Galanova, N., & Vedernikova, M. (2013). Testing Goodness of Fit of Parametric AFT and PH Models with Residuals. \emph{Communications in Statistics - Simulation and Computation}, 42(6), 1352–1367.

Wang H. J., and Wang, L. (2009). Locally weighted censored quantile regression. \emph{Journal of the American Statistical Association}, 104, 1117-1128. 

Martin, M. A., and Roberts, S. (2010). Jackknife-after-bootstrap regression influence diagnostics, \emph{Journal of Nonparametric Statistics}, 22, 257-269.

Nardi, A., and Schemper, M. (1999). New residuals for cox regression and their application to outlier screening, \emph{Biometrics}, 55, 523-529.
}
\note{
  We would like to thank Huxia Judy Wang and Lan Wang for the permission to use the function LCRQ().
}
\seealso{
    \code{\link{odc}}, \code{\link{plot}}, \code{\link{coef}}, \code{\link{show}}, \code{\link[quantreg:crq]{quantreg}}
}
\keyword{ package }
%end