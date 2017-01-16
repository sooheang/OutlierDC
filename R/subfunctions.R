.onAttach <- function(libname, pkgname) {
    ver <- read.dcf(file=system.file("DESCRIPTION", package=pkgname),
                      fields="Version")
    packageStartupMessage("")
    packageStartupMessage("Package ", pkgname, " (",ver,") loaded.")
}

##################################################

####################################################################
# A function bodc() to provide bootstrp standard errors for scoring algorithm
# for the scoring algorithm
bodc <- function(object, B = 500, q = 0.05, fence = c("UB", "LB", "both")){

  #####
	# Check class!!
	# pre-processing
	fence <- match.arg(fence)

	data = object@raw.data
	n = nrow(data)
	n.half = floor(n/2)
	reg = object@reg

	s = object@score
	s.ord = order(s)
	s.ord2 = order(s.ord)
	s.diff = diff(s[s.ord])

	outlier = rep(FALSE, n)

	# generate bootstrap samples
	boot.index = matrix(1:n, nrow = n, ncol = B)
	boot.index <- apply(boot.index, 2, 
	                    function(x, n) sample(x, size = n, replace = TRUE), 
	                    n = n)

	boot.dist = rep(NA, B)

	#jack.sample = boot.index[ ,!apply(boot.index, 2, function(x) any(x == i)), drop = FALSE]
	cat("Now bootstrapping.....\n")
	# NOTICE ME: have to reduce computational cost by using vectorized-computation
	for(j in 1:B){
		# Extend the JaB sample for other detection algorithms
		fit.tmp = odc(
		  formula = formula(object@formula), 
		  data = data[boot.index[,j],], 
		  alg = "score", 
		  reg = reg)

		if(all(fit.tmp@fitted.mat[,3] <= fit.tmp@fitted.mat[,5])){
			# to avoid quantile-crossing problem
			score.tmp = fit.tmp@score
			score.tmp <- score.tmp[!(score.tmp %in% c(NA, NaN, Inf, -Inf))]
			boot.dist[j] <- sd(score.tmp) 
		} 
	}

	#object@upper <- quantile(boot.dist, 0.01)
	#s.diff.out = which(s.diff > object@upper)
	object@ks <- quantile(boot.dist, q, na.rm = TRUE)
	
	# step 2. declare outliers
	if(fence %in% c("UB", "both")){
		#con = s.diff.out >= floor(quantile(1:n, probs = ifelse(fence == "UB", 1-alpha, 1 - (alpha / 2))))
		#out.UB = s.diff.out[con]
		out.UB = which(s.diff[n.half:n] > object@ks)
		if(length(out.UB) >= 1){
			out.UB <- names(out.UB[1])
			outlier[s.ord[which(s.ord == out.UB):n]] <- TRUE
		}
	} else if (fence %in% c("LB", "both")){
		#out.UB = s.diff.out[s.diff.out <= ceiling(quantile(1:n, probs = ifelse(fence == "LB", alpha, alpha / 2)))]
		out.LB = which(s.diff[1:n.half] > object@ks)
		if(length(out.LB) >= 1){
			out.LB <- names(out.LB[length(out.LB)])
			outlier[s.ord[1:which(s.ord == out.UB)]] <- TRUE
		}
	} 

	#declare outliers
	object@outliers <- outlier
	object@n.outliers <- sum(outlier)
	object@boot.dist <- boot.dist
	object@refined.data <- object@raw.data[!object@outliers,, drop = FALSE]
	object@boot.index <- boot.index
	cat("Done.\n")
	return(object)
}

####################################################################
# A function JaB() to provide the routine for a Jackknife-after-Bootstrap resampling 
# for the scoring algorithm
JaB <- function(object, B = 1000, alpha = 0.05, fence = c("UB", "LB", "both")){
	#object: OutlierDC object
	# i: a row mumber for the ith observation based on the rownames() 
	# plot: a logical paramter for the display of the JaB QQ plot
	# Jackknife-after-Bootstrap (JaB) with the scoring algorithm
	
	#####
	# Check class!!
	# pre-processing
	fence <- match.arg(fence)

	data = object@raw.data
	n = nrow(data)
	reg = object@reg
	
	# Determine which observation is an outlier 
	# by using (modified) Jackknife-after-Bootstrap (JaB) sampling
	
	# generate Jackknife after Bootstrap samples
	# i = 1; j = 1
	#scores.dist = list()
	scores.cri = matrix(0, nrow = n, ncol = 5)
	colnames(scores.cri) <- c("Score", "LB", "UB", "Criteria", "Quantile")
	#rownames(scores.cri) <- names(object@score)
	scores.cri <- as.data.frame(scores.cri)

	scores.cri[,1] <- object@score
	scores.ord = order(scores.cri[,1])
	scores.cri <- scores.cri[scores.ord,]

	# generate bootstrap samples
	#names.set = names(object@score)
	boot.index = matrix(1:n, nrow = n, ncol = B)
	boot.index <- apply(boot.index, 2, function(x, n) sample(x, size = n, replace = TRUE), n = n)

	score.boot = matrix(NA, nrow = n, ncol = B)
	#jack.sample = boot.index[ ,!apply(boot.index, 2, function(x) any(x == i)), drop = FALSE]

	# obtain the sampling distribution for the ith observation
	cat("Now bootstrapping.....\n")
	# NOTICE ME: have to reduce computational cost by using vectorized-computation
	for(j in 1:B){
		# Extend the JaB sample for other detection algorithms
		fit.tmp = odc(formula = formula(object@formula), data = data[boot.index[,j],], alg = "score", reg = reg)
		if(all(fit.tmp@fitted.mat[,3] <= fit.tmp@fitted.mat[,5])){
			# avoid quantile-crossing problem
			score.boot[,j] <- fit.tmp@score
			#score.ind <- c(score.ind, fit.tmp@score)
		} 
	}

	# estimated pth quantile function definition of the JaB distribution
	for(i in 1:n){
		ind.dist = score.boot[ ,apply(boot.index, 2, function(x,i.subset) any(x == i.subset), i.subset = i)]
		ind.dist <- as.vector(ind.dist)
		ind.dist <- ind.dist[!(ind.dist %in% c(Inf, -Inf, NA, NaN))]
		#ind.dist <- na.omit(ind.dist)
		if(fence == "UB"){
			scores.cri[i,3] <- quantile(ind.dist, 1 - alpha)
			scores.cri[i,4] <- scores.cri[i,1] >= scores.cri[i,3]
		} else if(fence == "LB"){
			scores.cri[i,2] <- quantile(ind.dist, alpha, na.rm = TRUE)
			scores.cri[i,4] <- object@score[i] <= scores.cri[i,2]
			#print(scores.cri[i,])
		} else if(fence == "both"){
			scores.cri[i,2] <- quantile(ind.dist, alpha / 2, na.rm = TRUE)
			scores.cri[i,3] <- quantile(ind.dist, 1 - alpha / 2, na.rm = TRUE)
			scores.cri[i,4] <- object@score[i] < scores.cri[i,2] | object@score[i] > scores.cri[i,3]
		}

		scores.cri[i,5] <- quantile(ind.dist, ifelse(i == 1, i / (n - 1), 
				ifelse(i == n, i / (n+1), (i - 1) / (n -1) )), na.rm = TRUE)
	}

	# set the subset of the alpha for the order statistic of the outlying score
	# alpha.upper = floor( quantile(1:n, 1-alpha) )
	# alpha.lower = ceiling( quantile(1:n, alpha) )
	
	#declare outliers
	scores.ord2 = order(scores.ord)
	object@score <- scores.cri[scores.ord2,1]
	object@lower <- scores.cri[scores.ord2,2]
	object@upper <- scores.cri[scores.ord2,3]
	object@outliers <- as.logical(scores.cri[scores.ord2,4])
	object@n.outliers <- as.integer(sum(object@outliers))
	object@boot.dist <- scores.cri[scores.ord2,5]
	object@refined.data <- object@raw.data[!object@outliers,, drop = FALSE]
	object@score.boot <- score.boot
	object@boot.index <- boot.index
	return(object)
}

##################################################
setGeneric("show")
setMethod("show", "OutlierDC", function(object){

		cat("\n     Outlier Detection for Censored Data\n\n")
		cat(" Call: ")
		print(object@call)
		
		cat("\n Algorithm: ")
		switch(object@alg, 
			score = cat("Scoring algorithm"),
			boxplot = cat("Boxplot algorithm"),
			residual = cat("Residual-based algorithm")
		)
		cat(paste(" (",object@alg,")", sep =""),"\n")

		cat(" Model: ")
		switch(object@reg, 
			#coxph = cat("Cox proportional hazard regresion"),
			PengHuang = cat("Censored quantile regression with the Peng and Huang estimator"),
			Powell = cat("Censored quantile regression with the Powell estimator"),
			Portnoy = cat("Censored quantile regression with the Portnoy estimator"),
			Wang = cat("Locally weighted censored quantile regression")
		)
		cat(paste(" (",object@reg,")", sep = ""),"\n")

		mf1 <- Formula::model.frame(object@formula, data = object@raw.data)
		resp <- Formula::model.response(mf1)
		times = resp[ ,1]
		delta = resp[ ,2]
		X = Formula::model.matrix(object@formula, data = object@raw.data)
		n = length(times)

		cat(" # of outliers detected: ", object@n.outliers, "\n")

		if(object@alg == "residual") {
			cat(" Value for cut-off kr: ",object@kr,"\n")
		} else if(object@alg == "boxplot") {
			cat(" Value for cut-off kb: ",object@kb,"\n")
		} else if(object@alg == "score") {
			cat(" Value for cut-off ks: ", object@ks,"\n")
		}

		
		if(object@alg == "residual"){
			print.data = cbind(times, delta, X[,-1], 
			                   residual = object@score, 
			                   cutoff = ifelse(object@reg == "coxph", 
			                                   object@cutoff, 
			                                   object@kr * object@cutoff)
			                   )

			order.row = order(times)
			print.data <- zapsmall(print.data, digits = 3)
			print.data <- as.data.frame(print.data[object@outliers, , drop = FALSE])
			order.row <- order(print.data$times, decreasing = TRUE)
			print.data <- cbind(print.data, Outlier = "*")
			#order <- order(object@score, decreasing = TRUE)
			#print.data <- zapsmall(print.data, digits = 3)
			#Signif <- ifelse(object@outliers, "*", "")
			#print.data <- cbind(print.data, Outlier = Signif)
			#print.data <- as.data.frame(print.data)
			#print.data <- print.data[order, , drop = FALSE]

			if(object@fence %in% c("both", "UB")){
				cat("\n Outliers detected:\n")
				print.head <- head(print.data)
				print(print.head)
				cat("\n", nrow(print.head) ,"of all", object@n.outliers, "outliers were displayed. \n")
			}

			if(object@fence %in% c("both", "LB")){
				cat("\n Outliers detected by lower fence:\n")
				#print(print.data[n:(n-5),], digits = 3)
				print(tail(print.data))
			}
		} else if(object@alg == "boxplot"){
			cat("\n Outliers detected:\n")
			if(object@fence == "UB"){
				print.data <- cbind(times, delta, X, UB = object@upper)
			}
			else if(object@fence == "LB"){
				print.data <- cbind(times, delta, X, LB = object@lower)				
			}
			else if(object@fence == "both"){
				print.data <- cbind(times, delta, X, LB = object@lower, UB = object@upper)				
			}

			order.row <- order(times)
			print.data <- zapsmall(print.data, digits = 3)
			print.data <- as.data.frame(print.data[object@outliers, , drop = FALSE])
			order.row <- order(print.data$times, decreasing = TRUE)
			print.data <- cbind(print.data, Outlier = "*")
			print.order <- print.data[order.row, ]
			print(print.order)
			cat("\n", nrow(print.order) ,"of all", object@n.outliers, "outliers were displayed. \n")
		} else if(object@alg ==  "score"){
			print.data <- cbind(T = times, C = delta, X[,-1, drop = FALSE], score = object@score)
			
			order <- order(object@score, decreasing = TRUE)
			print.data <- zapsmall(print.data, digits = 3)

			Signif <- ifelse(object@outliers, "*", "")
			print.data <- cbind(print.data, Outlier = Signif)
			print.data <- as.data.frame(print.data)
			print.data <- print.data[order, , drop = FALSE]

			if(object@fence %in% c("both", "UB")){
				cat("\n Top 6 outlying scores:\n")
				print(head(print.data))
			}
			
			if(object@fence %in% c("both", "LB")){
				cat("\n Bottom 6 outlying scores:\n")
				#print(print.data[n:(n-5),], digits = 3)
				print(tail(print.data))
			}						
		} 
	}
)

####################################################################
setGeneric("plot")
setMethod("plot", "OutlierDC", function(x, y = NA, ...){
		mf1 <- model.frame(x@formula, x@raw.data)
		resp <- model.response(mf1)
		Times <- resp[ ,1]
		status <- resp[ ,2]

		if(x@alg == "residual"){
			Residuals = x@score
			Fitted.values = x@fitted.mat[ ,4]
			limit.y <- max(Residuals, abs(x@kr) * x@cutoff)
			plot(Fitted.values, Residuals,pch = c(1,3)[status+1], ylim = c(-1 * limit.y, limit.y), ...)
			grid()
			points(Fitted.values[x@outliers], Residuals[x@outliers], pch = c(1,3)[status[x@outliers]+1], col = "blue")
			if(x@fence %in% c("both", "UB")) abline(h = x@kr * x@cutoff, col = "blue", lty = 2, lwd = 1.5)
			if(x@fence %in% c("both", "LB")) abline(h = -1 * x@kr * x@cutoff, col = "blue", lty = 2, lwd = 1.5)
			legend("bottomleft",c("Censored","Event"), cex=1, pch=c(1,3), bty = "n")
		}
		else if(x@alg == "score"){
			Scores = x@score
			if(as.logical(length(x@boot.dist))){
				# JaB Q-Q plot of outlying scores
				#plot(x@boot.dist, Scores, 
				#	main = "Empirical Q-Q plot of outlying scores",
				#	pch = c(1,3)[status+1],
				#	xlab = "Empirical score",
				#	ylab = "Outlying score"
				#)
					#xlim = range(x@JaB.dist),
					#ylim = range(x@JaB.dist)
				#)
				#points(x@boot.dist[x@outliers], Scores[x@outliers], col = "red", cex = 1, pch = c(1,3)[status[x@outliers] + 1])
				#legend("bottomright", c("Censored", "Event"), cex = 1, pch = c(1,3), bty = "n")
				#abline(c(0,1), col = "blue", lwd = 1.25, lty = 2)
			} else{
				# normal Q-Q plot of outlying scores
				tmp <- qqnorm(Scores, main = "Q-Q plot of outlying scores", pch = c(1,3)[status+1])
				qqline(Scores, col = "tomato", lwd = 1.5)
				#grid()
				if(!is.logical(x@upper)) abline(h = x@upper, col = "blue", lwd = 2, lty = 2)
				if(!is.logical(x@lower)) abline(h = x@lower, col = "blue", lwd = 2, lty = 2)

				points(tmp$x[x@outliers], tmp$y[x@outliers], col = "red", pch = c(1,3)[status[x@outliers]+1])
				legend("bottomright",c("Censored","Event"), cex=1, pch=c(1,3), bty = "n")
			}
		}
		else if(x@alg == "boxplot"){
			cov.x <- model.matrix(x@formula, data = x@raw.data)
			n <- ncol(cov.x) - 1
			for(i in 1:n){
				covariate <- cov.x[ ,i+1]
				order <- order(covariate)
				plot(covariate, Times, pch = c(1,3)[status+1], axes=F, ...)
				points(covariate[x@outliers], Times[x@outliers], pch =  c(1,3)[status[x@outliers]+1], col = "blue")
				
				if(x@fence %in% c("both", "UB")){
				lines(covariate[order], x@fitted.mat[order,5], lwd = 1.5,lty = 3)
				lines(covariate[order], x@upper[order], col = "blue", lwd = 1.5,lty = 2)
				}
				lines(covariate[order], x@fitted.mat[order,4], lwd = 1.5)
				if(x@fence %in% c("both", "LB")){
				lines(covariate[order], x@fitted.mat[order,3], lwd = 1.5,lty = 3)
				lines(covariate[order], x@lower[order], col = "blue", lwd = 1.5,lty = 2)
				}
				axis(1)
				axis(2, at= round(quantile(Times, probs = 0:5/5),1), labels= round(quantile(Times, probs = 0:5/5),1))
				box()
				rug(jitter(cov.x[order], amount = 0.01), ticksize = 0.01)
				legend("topright",c("Censored","Event"), cex=1, pch=c(1,3), bty = "n")
			}
		}
	}
)

###################################################################
setGeneric("coef")
setMethod("coef", "OutlierDC", function(object) round(object@coefficients,3) )

####################################################################
setGeneric("summary")
setMethod("summary","OutlierDC", function(object, taus = c(.1, .25, .5, .75, .9)){
		fit <- crq(object@formula, data = object@raw.data, method = object@alg)
		summary(fit, taus = taus)
	}
)

####################################################################
setGeneric("update")
setMethod("update","OutlierDC", function(object, alpha = 0.01, ks = NA, LB = NA){
# object: OutlierDC object
# UB, LB: sample quantiles
# This function is designed for the scoring algorithm

	if(length(object@score.boot) == 0){
		# update ks 
		Scores = object@score
		UB = ks
		if(is.na(UB) & is.na(LB)) stop("Please, update the object using the argument UB and LB")
		else if(!is.na(UB) & is.na(LB)){
			if(!is.logical(object@lower)) object@outliers <- Scores > UB | object@outliers
			else object@outliers <- Scores > UB
			object@ks <- UB
			object@upper <- UB
		}
		else if(is.na(UB) & !is.na(LB)){
			if(!is.logical(object@upper)) object@outliers <- Scores < LB | object@outliers	
			else object@outliers <-  Scores < LB
			
			object@lower <- LB
		}
		else{
			mf1 <- model.frame(object@formula, data = object@raw.data)
			resp <- model.response(mf1)
			times = resp[,1]
		
			object@outliers <- ifelse(times > object@fitted.mat[,4],
					Scores > UB, 
					Scores < LB)
			object@lower <- LB
			object@upper <- UB
			object@ks <- UB
		}
		object@n.outliers <- as.integer(sum(object@outliers))
		object@refined.data <- object@raw.data[!object@outliers,, drop = FALSE]
		return(object)
	} else{
		data = object@raw.data
		n = nrow(data)
		reg = object@reg	

		scores.cri = matrix(0, nrow = n, ncol = 5)
		colnames(scores.cri) <- c("Score", "LB", "UB", "Criteria", "Quantile")
		scores.cri <- as.data.frame(scores.cri)

		scores.cri[,1] <- object@score
		scores.ord = order(scores.cri[,1])
		scores.cri <- scores.cri[scores.ord,]

		score.boot = object@score.boot
		boot.index = object@boot.index

		for(i in 1:n){
			ind.dist = score.boot[ ,apply(boot.index, 2, function(x,i) any(x == i), i = i)]

			ind.dist <- as.vector(ind.dist)
			ind.dist <- ind.dist[!(ind.dist %in% c(Inf, -Inf, NA))]
		
			if(object@fence == "UB"){
				scores.cri[i,3] <- quantile(ind.dist, 1 - alpha, na.rm = TRUE)
				scores.cri[i,4] <- scores.cri[i,1] >= scores.cri[i,3]
				#print(scores.cri[i,])
			} else if(object@fence == "LB"){
				scores.cri[i,2] <- quantile(ind.dist, alpha, na.rm = TRUE)
				scores.cri[i,4] <- object@score[i] <= scores.cri[i,2]
			} else if(object@fence == "both"){
				scores.cri[i,2] <- quantile(ind.dist, alpha / 2, na.rm = TRUE)
				scores.cri[i,3] <- quantile(ind.dist, 1 - alpha / 2, na.rm = TRUE)
				scores.cri[i,4] <- object@score[i] < scores.cri[i,2] | object@score[i] > scores.cri[i,3]
			}
		}
	
		# adjust outliers based on the alpha
		#alpha.ind1 = as.integer( quantile(1:n, alpha) )
		#alpha.ind2 = as.integer( quantile(1:n, 1-alpha) )

		#if( as.logical(scores.cri[alpha.ind1, 4]) ) scores.cri[1:alpha.ind1, 4] <- TRUE
		#if( as.logical(scores.cri[alpha.ind2, 4]) ) scores.cri[alpha.ind2:n, 4] <- TRUE

		#declare outliers
		scores.ord2 = order(scores.ord)
		object@lower <- scores.cri[scores.ord2,2]
		object@upper <- scores.cri[scores.ord2,3]
		object@outliers <- as.logical(scores.cri[scores.ord2,4])
		object@n.outliers <- as.integer(sum(object@outliers))
		object@refined.data <- object@raw.data[!object@outliers,, drop = FALSE]
		return(object)
	}
	# END of update()
	# Oct 10, 2014 by Soo-Heang EO
} 
)
# End of subfunctions()
# Oct 10, 2014 by Soo-Heang EO

