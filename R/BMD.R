BMD <- function(object, bmr = 0.10, backg = 0, def = 'additional', eq = 'as.is', sigLev = 0.05, ci = 'CI', sav = FALSE){

#BMD <- function(object, bmr = 0.05, backg = 0, rtype = c('quantal', 'continuous'), def = c('additional', 'excess', 'relative', 'hybrid'), eq = FALSE, sigLev = 0.05, display = FALSE){

# eq  as.is,  default, or a list of Eqs
	
	## Calculating benchmark dose (BMD) and lower limit of benchmark dose (BMDL)
	## Selecting the Benchmark Response Level (BMR) (https://www.epa.gov/bmds/benchmark-dose-bmd-methods#bmr)
	## Dichotomous variables. Data on dichotomous variables are commonly presented as a fraction or percent of individuals that exhibit the given 
	## condition at a given dose or exposure level. Note that for modeling dichotomous data, one uses the exact counts. For such endpoints, 
	## normally we select probability density models like logistic, probit, Weibull, and so forth, whose predictions lie between zero and one for 
	## any possible dose, including zero.
	
	## Continuous variables. Data for continuous variables are often presented as means and SDs or SEs but may also be presented as a percent of 
	## control or some other standard. From a modeling standpoint, the most desirable form for such data is by individual. Unlike the usual 
	## situation for dichotomous variables, summarization of continuous variables results in a loss of information about the distribution 
	## of those variables. In addition, individual data is required when the intention is to use covariates in the analysis.
	## find non-acii in regexp mode [^\x00-\x7F]+
	
	###################################################
	getGap <- function(fct, x, paramHat, mse){
		jac <- jacobian(fct, x, paramHat)
		covPara <- mse * solve(t(jac) %*% jac)
		gap.PI <- sqrt(mse + diag(jac %*% covPara %*% t(jac))) # prediction intervals
		gap.CI <- sqrt(diag(jac %*% covPara %*% t(jac))) # confidence intervals
		list(PI = gap.PI, CI = gap.CI)
	}
	###################################################
	getBmd <- function(fct, x, paramHat, mse, bmrScaled, def, probT, ci, tol){
		# calculate bmd, bmdl, and bmdu
		rspnRange <- CEx(fct, paramHat, (c(0, Inf)))
		print(rspnRange)
		if(def == 'hybrid'){
			bmr0 <- min(rspnRange) + abs(rspnRange[1] - rspnRange[2]) * bmrScaled
		}else{
			bmr0 <- bmrScaled
		}

		if(bmr0 <= min(rspnRange)) stop('low bmr or high backg...')

		conc <- ECx(fct, paramHat, bmr0)
		xx <- c(conc, x)
		gap <- getGap(fct, xx, paramHat, mse)
		
		if(ci == 'PI'){
			bmrl <- bmr0 - probT * gap$PI[1]
			bmru <- bmr0 + probT * gap$PI[1]
		
		}else {
			bmrl <- bmr0 - probT * gap$CI[1]
			bmru <- bmr0 + probT * gap$CI[1]
		}
		bmrReal <- c(bmrl, bmr0, bmru)
		print(bmrReal)
		if(bmrReal[1] <= min(rspnRange)) bmrReal[1] <- min(rspnRange) * 1.0000001
		print(bmrReal)
		bmds <- ECx(fct, paramHat, bmrReal)
	}
	###################################################	

	tol <- 1E-18
	model <- object$eq
	param <- object$p
	crcInfo <- object$crcInfo
	size <- dim(crcInfo)
	n <- size[1]
	m <- length(param)
	x <- crcInfo[, 1]
	rspn <- crcInfo[, c(3 : (size[2] - 4))]
	rspn <- rowMeans(rspn)
	##
	rtype <- object$rtype
	res <- object$res
	resVar <- sum(res^2) / (n - m)
	
	###################################################		
	if(identical(rtype, 'quantal')){
		bmrScaled <- switch(def,
		"additional" = bmr * (1 - backg) + backg,
		"excess" = bmr + backg)
	}else if(identical(rtype, 'continuous')){
		rspnRange <- CEx(model, param, (c(0, Inf)))
		bmrScaled0 <- switch(def, 
		"relative" = bmr,
		# Kennyp, 2002. Critical Issues in Benchmark Calculations from Continuous Data. Crit. Rev. Toxicol. 32, 133-153.
		# Equation 9
		"hybrid" = sqrt(resVar) * (qnorm(1 - backg) - qnorm(1 - (backg + bmr))) / abs(rspnRange[1] - rspnRange[2])
		)
		bmrScaled <- as.numeric(bmrScaled0)
	}else{
		cat(rtype, ' dose-responses are not supported here')
	}
	###################################################	
	if(eq[1] == 'default'){
		if(identical(rtype, 'quantal')) {
			eqList <- c('Hill', 'Weibull', 'Logit')
			#eqList <- c('Hill', 'Weibull', 'Logit', 'Weibull_three', 'Hill_three', 'Logit_three')
		}else if (identical(rtype, 'continuous')) {
			eqList <- c('Weibull_four', 'Logit_four', 'Hill_four')
		}
	}else if (eq[1] == 'as.is'){
			eqList <- model
	}else {
		eqList <- eq
	}

	probT <- qt(1 - sigLev / 2, n - m)

	if(eq[1] == 'as.is'){
	
		fct <- eqList[1]
		paramHat <- param
		mse <- resVar
		print(fct)
		bmds <- getBmd(fct, x, paramHat, mse, bmrScaled, def, probT, ci, tol)
	
	}else{
		
		eqNum <- length(eqList)
		bmdMatrix <- matrix(0, eqNum, 3)
		
		for(j in seq_len(eqNum)) {
			tempfit <- tuneFit(x, rspn, eqList[j], bmrScaled)
			sta <- tempfit$sta
			staLen <- length(sta)
			mse <- sta[staLen - 4]^2
			paramHat <- sta[-c((staLen - 7) : staLen)]
			bmdMatrix[j, ] <- getBmd(eqList[j], x, paramHat, mse, bmrScaled, def, probT, ci, tol)
		}
		bmdMatrix[which(bmdMatrix == NaN)] <- tol
		bmds <- bmdMatrix[which(bmdMatrix[, 1] == min(bmdMatrix[, 1])), ]
		bmds <- t(bmds)
	}
	
	colnames(bmds) <- c('BMDL', 'BMD', 'BMDU')

	if (sav != FALSE){
		if(sav == TRUE) {
			svfile = paste("BMD_bmr_", bmr, "_", Sys.Date(), ".txt", sep = "")
			write.table(bmds, svfile, sep = "\t", quote = F)
		} else{
			write.table(bmds, sav, sep = "\t", quote = F)
		}
	}

	return(bmds)
	
}
