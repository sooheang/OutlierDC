##########################################################################
#
#    Outlier Detection for Censored Data 
#
#      by
#      Soo-Heang Eo, Ph.D. 
#      Deparment of Statistics 
#      Korea University, Seoul
#
#      May 2015
#
##########################################################################

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
