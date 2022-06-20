getCI <- function(object, effv, Scaled = TRUE, sigLev = 0.05, sav = FALSE){
	## object generated using curveFit
	## Scaled(TRUE/FALSE) was used only in continuous dose-response (rtype == 'continuous')
	##
	############################################
	getGap <- function(fct, x, paramHat, mse){
		jac <- jacobian(fct, x, paramHat)
		covPara <- mse * solve(t(jac) %*% jac)
		gap.PI <- sqrt(mse + diag(jac %*% covPara %*% t(jac))) # prediction intervals
		gap.CI <- sqrt(diag(jac %*% covPara %*% t(jac))) # confidence intervals
		list(PI = gap.PI, CI = gap.CI)
	}
	############################################	
	
	tol <- 1E-9
	maxtol <- 1E6
	rtype <- object$rtype
	fct <- object$eq
	param <- object$p
	rspnRange <- object$rspnRange
	crcInfo <- object$crcInfo
	size <- dim(crcInfo)
	n <- size[1]
	m <- length(param)
	x <- crcInfo[, 1]
	rspn <- crcInfo[, c(3 : (size[2] - 4))]
	mse <- sum(object$res^2) / (n - m)

	if(rtype == 'continuous' && Scaled == TRUE){
		#rspnRange <- CEx(fct, param, c(0, Inf))
		effv0 <- rspnRange[1] + (rspnRange[2] - rspnRange[1]) * effv
	}else{
		effv0 <- effv
		#rspnRange <- CEx(fct, param, c(0, Inf))
	}

	probT <- qt(1 - sigLev / 2, n - m)

	if(rtype != 'hormesis'){
		ecx <- ECx(fct, param, effv0, rtype, FALSE)
		if(is.list(ecx)) ecx <- ecx$ecx
		xmat <- matrix(0, length(effv0), 6)
		emat <- matrix(0, length(effv0), 5)
		
		for(j in seq(length(effv0))){
			xx <- c(ecx[j], x)
			gap <- getGap(fct, xx, param, mse)
			plow <- ifelse(effv0[j] - probT * gap$PI[1] <= rspnRange[1], rspnRange[1] + (rspnRange[2] - rspnRange[1]) / maxtol, effv0[j] - probT * gap$PI[1])
			pup <- effv0[j] + probT * gap$PI[1]
			clow <- ifelse(effv0[j] - probT * gap$CI[1] <= rspnRange[1], rspnRange[1] + (rspnRange[2] - rspnRange[1]) / maxtol, effv0[j] - probT * gap$CI[1])
			cup <- effv0[j] + probT * gap$CI[1]
			xmat[j, c(3 : 6)] <- ECx(fct, param, c(plow, pup, clow, cup), rtype, Scaled = FALSE)
			xmat[j, 1] <- effv0[j]
			xmat[j, 2] <- ecx[j]
			emat[j, 1] <- effv0[j]
			emat[j, -1] <- c(plow, pup, clow, cup)
		}
		
		effv0 <- as.numeric(format(effv0, digit = 3))
		
		if(rtype == 'quantal'){
			rownames(xmat) <- paste0('EC', effv0 * 100)
			rownames(emat) <- paste0('E', effv0 * 100)
		}else{
			rownames(xmat) <- paste0('EC', effv0)
			rownames(emat) <- paste0('E', effv0)
		}
		colnames(xmat) <- c('effv_abs', 'EC', 'PI.low', 'PI.up', 'CI.low', 'CI.up')
		colnames(emat) <- c('effv_abs', 'PI.low', 'PI.up', 'CI.low', 'CI.up')
	
	}else{
	## hormesis dose-response
		if(min(effv0) < 0){
		## with negative effect
			effv_neg <- effv0[effv0 < 0]
			effv_pos <- effv0[effv0 >= 0]
			effvNew <- c(effv_neg, effv_neg, effv_pos)
			len_neg <- length(effv_neg)
			len_pos <- length(effv_pos)
			len_new <- length(effvNew)
			
			emat <- matrix(0, len_new, 5)
			xmat <- matrix(0, len_new, 6)
			ecx <- nmECx(fct, param, effv0, object$minx, rspnRange[1])
			
			for(j in seq(len_new)){
				xx <- c(ecx[j], x)
				gap <- getGap(fct, xx, param, mse)
				plow <- ifelse(effvNew[j] - probT * gap$PI[1] <= object$miny, object$miny + tol, effvNew[j] - probT * gap$PI[1])
				pup <- effvNew[j] + probT * gap$PI[1]
				clow <- ifelse(effvNew[j] - probT * gap$CI[1] <= object$miny, object$miny + tol, effvNew[j] - probT * gap$CI[1])
				cup <- effvNew[j] + probT * gap$CI[1]
				emat[j, ] <- c(effvNew[j], plow, pup, clow, cup)
				xmat[j, 1] <- effvNew[j]
				xmat[j, 2] <- ecx[j]
			}
		
			for(i in seq(len_neg)){
				temp <- nmECx(fct, param, emat[i, -1], object$minx, rspnRange[1])
				xmat[i, -c(1, 2)] <- temp[1 : 4]	
			}
			
			for(i in seq(len_pos + len_neg)){
				temp <- nmECx(fct, param, emat[(i + len_neg), -1], object$minx, rspnRange[1])
				xmat[i + len_neg, -c(1, 2)] <- temp[(length(temp) - 3) : length(temp)]
			}
			
			left_name <- paste0('ECL', effv_neg * 100)
			right_name <- paste0('ECR', effv_neg * 100)
			if(len_pos > 0) pos_name <- paste('EC', effv_pos * 100, sep = '') else pos_name <- c()
			rownames(xmat) <- c(left_name, right_name, pos_name)
			
			left_name <- paste0('EL', effv_neg * 100)
			right_name <- paste0('ER', effv_neg * 100)
			if(len_pos > 0) pos_name <- paste('E', effv_pos * 100, sep = '') else pos_name <- c()
			rownames(emat) <- c(left_name, right_name, pos_name)
			colnames(xmat) <-  c('effv_abs', 'EC','PI.low', 'PI.up', 'CI.low', 'CI.up')
			colnames(emat) <-  c('effv_abs', 'PI.low', 'PI.up', 'CI.low', 'CI.up')
			
		}else{
		## with only positive effect
			len_pos <- length(effv0)
			len_neg <- 0
			emat <- matrix(0, len_pos, 5)
			xmat <- matrix(0, len_pos, 6)
			ecx <- nmECx(fct, param, effv0, object$minx, rspnRange[1])
			
			for(j in seq(len_pos)){
				xx <- c(ecx[j], x)
				gap <- getGap(fct, xx, param, mse)
				plow <- ifelse(effv0[j] - probT * gap$PI[1] <= object$miny, object$minytol + tol, effv0[j] - probT * gap$PI[1])
				pup <- effv0[j] + probT * gap$PI[1]
				clow <- ifelse(effv0[j] - probT * gap$CI[1] <= object$miny, object$miny + tol, effv0[j] - probT * gap$CI[1])
				cup <- effv0[j] + probT * gap$CI[1]
				emat[j, ] <- c(effv0[j], plow, pup, clow, cup)
				xmat[j, 1] <- effv0[j]
				xmat[j, 2] <- ecx[j]
			}
			
			for(i in seq(length(effv0))){
				temp <- nmECx(fct, param, emat[i, -1], object$minx, rspnRange[1])
				xmat[i, -c(1, 2)] <- temp[(length(temp) - 3) : length(temp)]			
			}
			
			rowName <- paste0('EC', effv0 * 100)
			rownames(xmat) <- rowName
			rowName <- paste0('E', effv0 * 100)
			rownames(emat) <- rowName
			colnames(xmat) <-  c('effv_abs', 'EC','PI.low', 'PI.up', 'CI.low', 'CI.up')
			colnames(emat) <-  c('effv_abs', 'PI.low', 'PI.up', 'CI.low', 'CI.up')
		}
	}
	Results <- list(xmat = xmat, emat = emat)

	if (sav != FALSE){
		if(sav == TRUE) {
			sav = paste("getCI_", Sys.Date(), ".txt", sep = "")
		}
		sink(sav)
		print(Results)
		sink()
	}

	return(Results)
}
