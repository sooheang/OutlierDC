#' Outlier detection for the analysis of censored data
#' 
#' Detecting outlying observations for censored data, especially with an
#' application of lifetime studies.
#' 
#' The \code{odc} function conducts three outlier detection algorithms on the
#' basis of censored quantile regression. Three outlier detection algorithms
#' were implemented: residual-based, boxplot, and scoring algorithms. The
#' residual-based algorithm detects outlying observations using constant scale
#' estimates; however, it does not account for the heterogeneity of
#' variability. When the data is extremely heterogeneous, the boxplot algorithm
#' with censored quantile regression is more effective. The residual-based and
#' boxplot algorithms produce cut-offs to determine whether observations are
#' outliers. In contrast, the scoring algorithm provides the outlying magnitude
#' or deviation of each point from the distribution of observations. Outlier
#' detection is achieved by visualising the scores.
#' 
#' @param formula type of \code{Formula} object with a \code{survival} object
#' on the left-hand side of the ~ operator and covariate terms on the
#' right-hand side. The survival object with survival time and its censoring
#' status is constructed by the \code{\link[survival]{Surv}} function in
#' \code{survival} package.
#' @param data data frame with variables used in the \code{formula}. It needs
#' at least three variables, including survival time, censoring status, and
#' covariates.
#' @param alg type of an outlier detection algorithm. Three algorithms are
#' provided. The options \code{"score"}, \code{"boxplot"}, and
#' \code{"residual"} implement the scoring, boxplot, and residual-based
#' algorithms, respectively. The default algorithm is the scoring algorithm
#' with an option \code{"score"}.
#' 
#' @param reg type of a regression method used as a basis for outlier detection
#' algorithms. The options \code{"Wang"}, \code{"Portnoy"}, and
#' \code{"PengHuang"} conduct Wang and Wang's, Portnoy's, and Peng and Huang's
#' censored quantile regression approaches, respectively. The default is a
#' local CQR with an option \code{"Wang"}.
#' @param fence type of an outlying fence. Three options are provided. The
#' option \code{"UB"} provides the upper fence to find an observation who
#' "lived far too long". The option \code{"LB"} provides the lower fence to
#' find an observation who "died far too early". The option \code{"both"}
#' provides the upper and lower fences, simultaneously.
#' @param kr numeric value to control the tightness of cut-offs for the
#' residual algorithm with a default value of 3.
#' @param kb numeric value to control the tightness of cut-offs for the boxplot
#' algorithm with a default value of 1.5.
#' @param h bandwidth for locally weighted censored quantile regression with a
#' default value of 0.05.
#' @return an object of the S4 class "OutlierDC" with the following slots:\cr
#' call: evaluated function call\cr formula: formula to be used\cr raw.data:
#' data to be used for model fitting\cr refined.data: the data set after
#' removing outliers\cr refined.data: the data set containing outliers\cr
#' coefficients: the estimated censored quantile regression coefficient matrix
#' consisting of 10th, 25th, 50th, 75th, and 90th quantiles\cr fitted.mat: the
#' censored quantile regression fitted value matrix consisting of 10th, 25th,
#' 50th, 75th, and 90th quantiles\cr score: outlying scores (scoring algorithm)
#' or residuals (residual-based algorithm)\cr cutoff: estimated scale parameter
#' for the residual-based algorithm\cr lower: lower fence vector used for the
#' boxplot and scoring algorithms\cr upper: upper fence vector used for the
#' boxplot and scoring algorithms\cr outliers: logical vector to determine
#' which observations are outliers\cr n.outliers: number of outliers
#' detected\cr method: outlier detection method to be used\cr rq.model:
#' censored quantile regression to be used\cr kr: a value to be used for the
#' tightness of cut-offs in the residual algorithm \cr kb: a value to be used
#' for the tightness of cut-offs in the boxplot algorithm \cr ks: a value to be
#' used for the tightness of upper fence cut-offs used for the scoring
#' algorithm with the \code{\link{update}} function\cr fence: type of fence to
#' be used in the model fitting \cr alpha: numeric value for the significance
#' level \cr boot.dist: empirical quantiles by Jackknife-after-Bootstrapping
#' %B: a value to be used for the scoring algorithm in order to determine the
#' empirical distribution of each observation.  %LB: lower fence used for the
#' scoring algorithm with the \code{\link{update}} function\cr
#' @seealso \code{\link{OutlierDC-package}} \cr \code{\link{coef}},
#' \code{\link{plot}}, \code{\link{show}}, \code{\link{update}}
#' @source Eo, S-H, Hong, S-M, and Cho, H. (2014+). Identification of outlying
#' observations with quantile regression for censored data, \emph{Submitted}.
#' 
#' Martin, M. A., and Roberts, S. (2010). Jackknife-after-bootstrap regression
#' influence diagnostics, \emph{Journal of Nonparametric Statistics}, 22,
#' 257-269.
#' 
#' Nardi, A., and Schemper, M. (1999). New residuals for Cox regression and
#' their application to outlier screening, \emph{Biometrics}, 55, 523-529.
#' 
#' Wang H. J., and Wang, L. (2009). Locally weighted censored quantile
#' regression. \emph{Journal of the American Statistical Association}, 104,
#' 1117-1128.
#' @keywords odc
#' @examples
#' 
#'   \dontrun{
#' 
#' library(OutlierDC)
#' data(ebd)
#' str(ebd)
#' 
#' ####
#' # outlier detection using the scoring algorithm
#' fit = odc(Surv(log(time), status) ~ meta, data = ebd)
#' fit
#' 
#' # A threshold is added by k_s to this plot using the updata() function
#' fit1 = update(fit, ks = 4)
#' fit1
#' plot(fit1)
#' 
#' # A threshold can be determined by using the empirical distribution of standard deviation for the outlying scores
#' fit2 = bodc(fit, B = 500)
#' fit2
#' 
#' ####
#' # outlier detection using the residual-based algorithm
#' fit3 = odc(Surv(log(time), status) ~ meta, data = ebd, alg = "residual")
#' fit3
#' plot(fit3, main = "Residual-based algorithm")
#' 
#' 
#' fit3@outlier.data
#' 
#' ####
#' # outlier detection using the boxplot-based algorithm
#' fit4 = odc(Surv(log(time), status) ~ meta, data = ebd, alg = "boxplot")
#' fit4
#' plot(fit4, main = "Boxplot-based algorithm", xlab = "Number of metastatic lymph nodes",, ylab = "Log of survival times")
#' }
#' 
odc <- function(formula, data, 
	alg = c("score", "boxplot","residual"), 
	reg = c("Wang", "PengHuang", "Portnoy"), 
	fence = c("UB", "LB", "both"),
	kr = 3, kb = 1.5, h = .05){

  ##########
  #Preparation
  call <- match.call()
  alg <- match.arg(alg)
  reg <- match.arg(reg)
  fence <- match.arg(fence)

  if(!Formula::is.Formula(formula)) formula <- Formula::Formula(formula)

	if(is.data.frame(data)) data <- as.data.frame(data)

	mf1 = stats::model.frame(formula, data = data)
	X.mat = stats::model.matrix(formula, data = mf1)
	resp = stats::model.response(mf1)

	y = resp[,1]
	status = resp[,2]
  n = nrow(data)
  rownames(data) <- 1:n

	# NOTICE ME: need to add parallel computing to reduce computing cost
	# cat("Please wait... \n")
	betas = matrix(NA, nrow = ncol(X.mat), ncol = 7)

	result = methods::new("OutlierDC")
	result@call <- call
	result@formula <- formula

	#fit = coxph(formula(formula), method = "breslow", model = TRUE, data = data)
	if(reg == "Wang"){
    betas[,1] <- LCRQ(y = y, x = X.mat[ ,-1, drop = FALSE], delta = status, tau = .05, h= h)
		betas[,2] <- LCRQ(y = y, x = X.mat[ ,-1, drop = FALSE], delta = status, tau = .10, h= h)
  	betas[,3] <- LCRQ(y = y, x = X.mat[ ,-1, drop = FALSE], delta = status, tau = .25, h= h)
  	betas[,4] <- LCRQ(y = y, x = X.mat[ ,-1, drop = FALSE], delta = status, tau = .50, h= h)
  	betas[,5] <- LCRQ(y = y, x = X.mat[ ,-1, drop = FALSE], delta = status, tau = .75, h= h)
  	betas[,6] <- LCRQ(y = y, x = X.mat[ ,-1, drop = FALSE], delta = status, tau = .90, h= h)
  	betas[,7] <- LCRQ(y = y, x = X.mat[ ,-1, drop = FALSE], delta = status, tau = .95, h= h)
	} else if (reg %in% c("Portnoy", "PengHuang")) {
		fit = crq(formula, data = data, method = reg)
  	betas[,1] <- coef(fit, taus = .05)
  	betas[,2] <- coef(fit, taus = .10)
  	betas[,3] <- coef(fit, taus = .25)
		betas[,4] <- coef(fit, taus = .50)
  	betas[,5] <- coef(fit, taus = .75)
  	betas[,6] <- coef(fit, taus = .90)
  	betas[,7] <- coef(fit, taus = .95)
	} 

	fit.q05 = (X.mat %*% betas[,1])[,1]
	fit.q10 = (X.mat %*% betas[,2])[,1]
	fit.q25 = (X.mat %*% betas[,3])[,1]
	fit.q50 = (X.mat %*% betas[,4])[,1]
	fit.q75 = (X.mat %*% betas[,5])[,1]
	fit.q90 = (X.mat %*% betas[,6])[,1]
	fit.q95 = (X.mat %*% betas[,7])[,1]
	
	fitted.mat = cbind(fit.q05, fit.q10, fit.q25, fit.q50, fit.q75, fit.q90, fit.q95)
	result@fitted.mat <- fitted.mat

	############
	## Fit outlier detection algorithms
	if(alg == "score"){
	    
    # Step 1. calculate an outlying score
		score = ifelse(y >= fit.q50, 
  		(y - fit.q50) / (fit.q75 - fit.q50), 
			(y - fit.q50) / (fit.q50 - fit.q25))
		
		outlier = rep(FALSE, n)

		#s.ord = order(score)
		#s.ord2 = order(s.ord)
		#s.diff = diff(score[s.ord])
		#s.diff.out = which(s.diff > sd(score) * ks)
		# step 2. declare outliers
		#if(length(s.diff.out) >= 1){
		#	if(fence %in% c("UB", "both")){
		#		out.UB = s.diff.out[s.diff.out >= floor(quantile(1:n, probs = ifelse(fence == "UB", 1-alpha, 1 - (alpha / 2))))]
		#		if(length(out.UB) >= 1){
		#			out.UB <- out.UB[1]
		#			outlier[s.ord[(out.UB+1):n]] <- TRUE
		#		}

		#	} else if (fence %in% c("LB", "both")){
		#		out.UB = s.diff.out[s.diff.out <= ceiling(quantile(1:n, probs = ifelse(fence == "LB", alpha, alpha / 2)))]
		#		if(length(out.LB) >= 1){
		#			out.LB <- out.LB[length(out.LB)]
		#			outlier[s.ord[1:(out.LB)]] <- TRUE
		#		}
		#	} 
		#}

		n.outliers = sum(outlier)
		result@score <- score


		#result@score <- score
		#n.outliers <- as.integer(0)
	} else if(alg == "boxplot"){
  	# calculate fences
  	iqr <- fit.q75 -fit.q25
  	upper.fence <- fit.q75 + (kb * iqr)
  	lower.fence <- fit.q25 - (kb * iqr)

  	if(fence == "both"){
    	outlier <- ifelse(y > fit.q50,
    				y > upper.fence, y < lower.fence)
    	#outlier <- ifelse(is.na(outlier), FALSE, outlier)
  	} else if(fence == "LB"){
  		# calculate only lower fence
    	outlier <- y < lower.fence
    	#outlier <- ifelse(is.na(outlier), FALSE, outlier)
  	} else if(fence == "UB"){
  		# calculate only upper fence
    	outlier <- y > upper.fence
  	}

  	outlier <- ifelse(is.na(outlier), FALSE, outlier)
		result@lower <- lower.fence
		result@upper <- upper.fence
		n.outliers <- sum(outlier)

  } else if(alg == "residual"){

  	res.fit = y - fit.q50
  	#score <- score
  	cutoff = stats::median(abs(res.fit) / stats::qnorm(0.75), na.rm = TRUE)

  	if (fence == "both") {
  		outlier <- abs(res.fit) > (kr * cutoff)
  	} else if (fence == "UB") {
  		outlier <- res.fit > (kr * cutoff)
  	} else if (fence == "LB") {
  		outlier <- res.fit < -1 * (kr * cutoff)
  	}

  	#if (resid == "ndr"){
  		# Normal-deviate residuals
  		# Use suvefit.coxph() to obtain the estimated survival function of the fitted model, S(t,z)
  		# Ni = Phi^{-1}(S(t,z)) where Phi^{-1} is qnorm(S(t.z), mean = 0, sd = 1)

    	#	res.fit = survival::survexp(formula(formula), ratetable = fit, data = data, method = "individual.s")
    	#	res.fit <- res.fit - 1e-06
    	#	res.fit <- qnorm(res.fit)

    #	if (fence == "both") {
    #		cutoff = qnorm(1 - alpha / 2)
    #		outlier <- abs(res.fit) > cutoff
    #	} else if (fence == "UB") {
    #		cutoff = qnorm(1 - alpha)
    #		outlier <- res.fit > cutoff
    #	} else if (fence == "LB") {
    #		cutoff = qnorm(alpha)
    #		outlier <- res.fit < cutoff
    #	}
  	outlier <- ifelse(is.na(outlier), FALSE, outlier)
  	result@score <- res.fit # absolute value of residuals
  	result@cutoff <- cutoff
  	n.outliers <- sum(outlier)
  }

	## Declare outliers using D statistics
	## manipulate results object
	rownames(betas) <- colnames(X.mat)
	colnames(betas) <- c("q05", "q10", "q25", "q50", "q75", "q90", "q95")
	betas <- as.data.frame(betas)   
	
	rownames(data) <- names(y)
	result@raw.data <- data
	result@coefficients <- betas
	#result@resid <- resid
	result@alg <- alg
	result@reg <- reg
	result@outliers <- outlier	
	result@n.outliers <- n.outliers
	result@refined.data <- data[!outlier,, drop = FALSE]
	result@outlier.data <- data[outlier,, drop = FALSE]
	result@kr <- kr
	result@kb <- kb
	result@fence <- fence	
	#cat("Done. \n")
	return(result)
}
#END################################################################s
