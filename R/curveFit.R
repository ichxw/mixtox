curveFit <- function(x, rspn, eq, param, effv, rtype = 'quantal', sigLev = 0.05, sav = FALSE, ...){
	# NLS curve fitting for monotonic and non-monotonic equations
	# x a vector of treatment concentration
	# rspn a vector or matrix
	# for non-monotonic curve fitting, Brain_Consens, BCV, Hill_five, and Biphasic are highly recommended.
	# rtype response type: quantal, continuous, hormesis

	
	#############################################################
	## source('ECx.R')
	
	if (missing(x) || missing(rspn) || missing(eq) || missing(param)) stop('argument missing')
	#print('Please set the right rtype: quantal, continuous, hormesis')
	cat('\nProcessing', rtype, 'quantal dose-response data\n')
	n <- length(x) # the number of concentrations
	mode(param) <- "numeric"
	m <- length(param) # the number of parameters
	
	if (is.vector(rspn)){
		if (n != length(rspn)) stop("concentration and response should be in the same length")
		rspn <- as.matrix(rspn)
		y <- rspn
	}else if (is.matrix(rspn)){
		size <- dim(rspn)
		y <- rowMeans(rspn)
		if(n != size[1]) stop("concentration and response should be paired")
	}
	
	# non-monotonic or monotonic
	if(eq == 'Brain_Consens' || eq == 'BCV' || eq == 'Cedergreen' || eq == 'Beckon' || eq == 'Biphasic' || eq == 'Hill_five') Hormesis <- TRUE else Hormesis <- FALSE
	
	# define equation expression
	fun <- switch(eq,
		# For Hill equation: Alpha = EC50; Beta = m(Hill coefficient); Gamma = Top; Delta = Bottom
		Hill = 'y ~ 1 / (1 + (Alpha / x)^Beta)',
		# Howard GJ, Webster TF. 2009. Generalized concentration addition: A method for examining mixtures containing partial agonists. J. Theor. Biol. 259:469~477
		# Hill function with slope parameter 1. Alpha is EC50 here.
		Hill_two = 'y ~ Beta * x / (Alpha + x)',
		Hill_three = 'y ~ Gamma /(1 + (Alpha / x)^Beta)',
		Hill_four = 'y ~ Delta + (Gamma - Delta) / (1 + (Alpha / x)^Beta)',
		Hill_five = 'y ~ 1 - (1 + (Gamma - 1) / (1 + (Alpha / x)^Beta)) * (1 - 1 / (1 + (Delta / x)^Epsilon))',
		# Hill_nine = 'y ~ (Gamma / (1 + (Alpha / x)^Beta)) * (Gamma_one / (1 + (Delta / x)^Epsilon)) * (Gamma_two / (1 + (Alpha_two / x)^Beta_two))',
		Weibull = 'y ~ 1 - exp(-exp(Alpha + Beta * log10(x)))',
		Weibull_three = 'y ~ Gamma * (1 - exp(-exp(Alpha + Beta * log10(x))))',
		Weibull_four = 'y ~ Gamma + (Delta - Gamma) * exp(-exp(Alpha + Beta * log10(x)))',		
		Logit = 'y ~ 1/(1 + exp((-Alpha) - Beta * log10(x)))',
		Logit_three = 'y ~ Gamma / (1 + exp((-Alpha) - Beta * log10(x)))',
		Logit_four = 'y ~ Delta + (Gamma - Delta) / (1 + exp((-Alpha) - Beta * log10(x)))',
		BCW = 'y ~ 1 - exp(-exp(Alpha + Beta * ((x^Gamma - 1) / Gamma)))',
		BCL = 'y ~ (1 + exp(-Alpha - Beta *((x^Gamma - 1) / Gamma)))^(-1)',
		GL = 'y ~ 1 / (1 + exp(-Alpha - Beta * log10(x)))^Gamma',
		# An equation to describe dose responses where there isstimulatin of growth at low doses. 1989. Weed Research.
		Brain_Consens = 'y ~ 1 - (1 + Alpha * x) / (1 + exp(Beta * Gamma) * x^Beta)',
		# Vanewijk, P. H. and Hoekstra, J.A. Calculation of the EC50 and its confidence interval when subtoxic stimulus is present. 1993, Ecotoxicol. Environ. Saf.
		BCV = 'y ~ 1 - Alpha * (1 + Beta * x) / (1 + (1 + 2 * Beta * Gamma) * (x / Gamma)^Delta)',
		# Cedergreen, N., Ritz, C., Streibig, J.C., 2005. Improved empirical models describing hormesis. Environ. Toxicol. Chem. 24, 3166~3172
		Cedergreen = 'y ~ 1 - (1 + Alpha * exp(-1 / (x^Beta))) / (1 + exp(Gamma * (log(x) - log(Delta))))',
		# Beckon, W. et.al. 2008. A general approach to modeling biphasic relationships. Environ. Sci. Technol. 42, 1308~1314.
		Beckon = 'y ~ (Alpha + (1 - (Alpha) / (1 + (Beta / x)^Gamma))) / (1 + (x / Delta)^Epsilon)',
		# Zhu X-W, et.al . 2013. Modeling non-monotonic dose-response relationships: Model evaluation and hormetic quantities exploration. Ecotoxicol. Environ. Saf. 89:130~136;
		Biphasic = 'y ~ Alpha - Alpha / (1 + 10^((x - Beta) * Gamma)) + (1 - Alpha) / (1 + 10^((Delta - x) * Epsilon))'
	)
	
	## checking minpack.lm package, use the nlsLM or built-in nls for curve fitting
	#if(require(minpack.lm)){
	dframe <- data.frame(x, y)
	
	if(requireNamespace("minpack.lm", quietly = TRUE)){
		# print("use the minpack.lm package")
		
		if(eq == "Weibull" || eq == "Logit" || eq == "Hill" || eq == "Hill_two"){
			fit <- minpack.lm::nlsLM(fun, data = dframe, start = list(Alpha = param[1], Beta = param[2]), control = nls.lm.control(maxiter = 1000), ...)
			#fit <- minpack.lm::nlsLM(fun, data = dframe, start = list(Alpha = param[1], Beta = param[2]))
		}else if(eq == "BCW" || eq == "BCL" || eq == "GL" || eq == 'Brain_Consens' || eq == "Hill_three" || eq == 'Weibull_three' || eq == 'Logit_three'){
			fit <- minpack.lm::nlsLM(fun, data = dframe, start = list(Alpha = param[1], Beta = param[2], Gamma = param[3]), control = nls.lm.control(maxiter = 1000), ...)
		
		}else if(eq == 'BCV'|| eq == 'Cedergreen' || eq == "Hill_four" || eq == 'Weibull_four' || eq == 'Logit_four'){
			#fit <- minpack.lm::nlsLM(fun, data = dframe, start = list(Alpha = param[1], Beta = param[2], Gamma = param[3], Delta = param[4]), ...)
			fit <- minpack.lm::nlsLM(fun, data = dframe, start = list(Alpha = param[1], Beta = param[2], Gamma = param[3], Delta = param[4]), control = nls.lm.control(maxiter = 1000), ...)
		}else if(eq == 'Beckon' || eq == 'Biphasic' || eq == 'Hill_five'){
			fit <- minpack.lm::nlsLM(fun, data = dframe, start = list(Alpha = param[1], Beta = param[2], Gamma = param[3], Delta = param[4], Epsilon = param[5]), control = nls.lm.control(maxiter = 1000), ...)
			#fit <- minpack.lm::nlsLM(fun, data = dframe, start = list(Alpha = param[1], Beta = param[2], Gamma = param[3], Delta = param[4], Epsilon = param[5]))
		}else{
			stop('input the right equation name')
		}
		#detach(package: minpack.lm)
	
	}else {
		warning('please install package minpack.lm')
		if(eq == "Weibull" || eq == "Logit" || eq == "Hill" || eq == "Hill_two"){
			fit <- nls(fun, data = dframe, start = list(Alpha = param[1], Beta = param[2]), control = nls.lm.control(maxiter = 1000))
		
		}else if(eq == "BCW" || eq == "BCL" || eq == "GL" || eq == 'Brain_Consens' || eq == "Hill_three" || eq == 'Weibull_three' || eq == 'Logit_three'){
			fit <- nls(fun, data = dframe, start = list(Alpha = param[1], Beta = param[2], Gamma = param[3]), control = nls.lm.control(maxiter = 1000))
		
		}else if(eq == 'BCV'|| eq == 'Cedergreen' || eq == "Hill_four" || eq == 'Weibull_four' || eq == 'Logit_four'){
			#fit <- nls(fun, data = dframe, start = list(Alpha = param[1], Beta = param[2], Gamma = param[3], Delta = param[4]), control = nls.lm.control(maxiter = 1000))
			fit <- nls(fun, data = dframe, start = list(Alpha = param[1], Beta = param[2], Gamma = param[3], Delta = param[4]), control = nls.lm.control(maxiter = 1000))
		}else if(eq == 'Beckon' || eq == 'Biphasic' || eq == 'Hill_five'){
			fit <- nls(fun, data = dframe, start = list(Alpha = param[1], Beta = param[2], Gamma = param[3], Delta = param[4], Epsilon = param[5]), control = nls.lm.control(maxiter = 1000))
			#fit <- nls(fun, data = dframe, start = list(Alpha = param[1], Beta = param[2], Gamma = param[3], Delta = param[4], Epsilon = param[5]))
		}else{
			stop('input the right equation name')
		}
	}
	
	fitInfo <- summary(fit) # fitting information
	yhat <- predict(fit, x) # y prediction
	res <- as.vector(residuals(fit))
	sst <- sum((y - mean(y))^2) # total sum of squares
	sse <- sum(res^2) # sum of squared errors
	r2 <- 1 - sse / sst # coefficient of determination
	adjr2 <- 1 - sse * (n - 1) / (sst * (n - m)) # adjusted coefficient of determination
	rmse <- sqrt(sse / (n - m)) # root-mean-square error
	mae <- sum(abs(res)) / n # mean absolute error
	#Spiess A-N, Neumeyer N. 2010. An evaluation of R2 as an inadequate measure for nonlinear models in pharmacological and biochemical research: A Monte Carlo approach. BMC Pharmacol. 10: 11.
	lnL <- 0.5 * (-n * (log(2 * pi) + 1 - log(n) + log(sse)))
	aic <- 2 * m - 2 * lnL # Akaike information criterion 
	aicc <- aic + 2 * m * (m + 1) / (n - m - 1)
	bic <- m * log(n) - 2 * lnL # Bayesian information criterion
	sta <- t(c(r2, adjr2, mae, rmse, aic, aicc, bic))
	colnames(sta) <- c('r2', 'adjr2', 'MAE', 'RMSE', 'AIC', 'AICc', 'BIC')
	
	paramHat <- coef(fit)	
	jac <- jacobian(eq, x, paramHat) # jacobian matrix calculation
	probT <- qt(1 - sigLev / 2, n - m) # the student t distribution
	mse <- rmse^2  # squared residual standard error
	covPara <- mse * solve(t(jac) %*% jac)  # covariance matrix of the parameter estimates
	
	gap.PI <- sqrt(mse + diag(jac %*% covPara %*% t(jac))) # prediction intervals
	gap.CI <- sqrt(diag(jac %*% covPara %*% t(jac))) # confidence intervals
	
	PI.up <- yhat + probT * gap.PI # PI upper bound
	PI.low <- yhat - probT * gap.PI # PI lower bound
	CI.up <- yhat + probT * gap.CI # CI upper bound
	CI.low <- yhat - probT * gap.CI # CI lower bound
	crcInfo <- cbind(x, yhat, rspn, PI.low, PI.up, CI.low, CI.up)
	
	# compute highest stimulation (minimum effect) of the J-shaped curve and associated concentration. 
	# Brain_Consens, BCV, Cedergreen, Beckon, Biphasic
	if(Hormesis == TRUE){
		rtype <- 'hormesis'
		if(eq == 'Brain_Consens') Alpha = paramHat[1]; Beta = paramHat[2]; Gamma = paramHat[3]
		if(eq == 'BCV' || eq == 'Cedergreen') Alpha <- paramHat[1]; Beta <- paramHat[2]; Gamma <- paramHat[3]; Delta <- paramHat[4]
		if(eq == 'Beckon' || eq == 'Biphasic' || eq == 'Hill_five') Alpha <- paramHat[1]; Beta <- paramHat[2]; Gamma <- paramHat[3]; Delta <- paramHat[4]; Epsilon <- paramHat[5]
	
		if (eq == 'Brain_Consens') f <- function(x) 1 - (1 + Alpha * x) / (1 + exp(Beta * Gamma) * x^Beta)
		if(eq == 'BCV') f <- function(x) 1 - Alpha * (1 + Beta * x) / (1 + (1 + 2 * Beta * Gamma) * (x / Gamma)^Delta)
		if(eq == 'Cedergreen') f <- function(x) 1 - (1 + Alpha * exp(-1 / (x^Beta))) / (1 + exp(Gamma * (log(x) - log(Delta))))
		if(eq == 'Beckon') f <- function(x) (Alpha + (1 - (Alpha) / (1 + (Beta / x)^Gamma))) / (1 + (x / Delta)^Epsilon)
		if(eq == 'Biphasic') f <- function(x) Alpha - Alpha / (1 + 10^((x - Beta) * Gamma)) + (1 - Alpha) / (1 + 10^((Delta - x) * Epsilon))
		if(eq == 'Hill_five') f <- function(x) 1 - (1 + (Gamma - 1) / (1 + (Alpha / x)^Beta)) * (1 - 1 / (1 + (Delta / x)^Epsilon))
		
		# intervals for finding the minimum
		intv <- c(x[1], x[length(x) - 1])

		minxy <- tryCatch({
			minxy <- optimize(f, intv)
		}, warning = function(w){
			message("Input an optimal intv")
		}, finally = {
			minxy <- list(minimum = NULL, objective = NULL)
		})
		minx <- minxy$minimum
		miny <- minxy$objective
	}

	## checking argument
	if(!missing(effv)){
		if(rtype == 'quantal' || rtype == 'continuous'){
			## effect concentration and confidence intervals 
			ecx <- ECx(eq, paramHat, effv, rtype = rtype)
		}else if (rtype == 'hormesis'){
			effv <- sort(effv)
			ecx <- nmECx(eq, paramHat, effv, minx)
		}
	}else{
		ecx <- NULL
	}
	
	#rspnRange <- CEx(eq, paramHat, c(0, Inf))
	rspnRange <- CEx(eq, paramHat, c(0, 1e30))
	
	if(Hormesis == FALSE){
		if(is.list(ecx)){
			Results <- list(fitInfo = fitInfo, eq = eq, p = paramHat, res = res, sta = sta, crcInfo = crcInfo, effvAbs = ecx$effvAbs, ecx = ecx$ecx, rtype = rtype, rspnRange = rspnRange)
		}else{
			Results <- list(fitInfo = fitInfo, eq = eq, p = paramHat, res = res, sta = sta, crcInfo = crcInfo, ecx = ecx, rtype = rtype, rspnRange = rspnRange)
		}
	}else{
		Results <- list(fitInfo = fitInfo, eq = eq, p = paramHat, res = res, sta = sta, minx = minx, miny = miny, crcInfo = crcInfo, ecx = ecx, rtype = rtype, rspnRange = rspnRange)
	}

	if (sav != FALSE){
		if(sav == TRUE) {
			sav = paste("curveFit_", eq, "_",Sys.Date(), ".txt", sep = "")
		}
		sink(sav)
		print(Results)
		sink()
	}

	return(Results)
}
