setOldClass("Formula")
setOldClass("Surv")
setClass(Class = "OutlierDC", representation(
							call 		= "language",
							formula 	= "Formula",
							raw.data 	= "data.frame",
							refined.data = "data.frame",
							outlier.data = "data.frame",
							coefficients ="data.frame",
							fitted.mat 	= "matrix",
							score 		= "vector",
							score.boot	= "matrix",
							boot.index 	= "matrix",
							boot.dist	= "numeric",
							cutoff 		= "vector",
							lower 		= "vector",
							upper 		= "vector",
							outliers 	= "vector",
							n.outliers 	= "integer",
							alg 		= "character",
							reg 		= "character",
							kr 			= "numeric",
							kb 			= "numeric",
							ks 			= "numeric",
							fence 		= "character",
							alpha 		= "numeric")
)
