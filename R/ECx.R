ECx <- function(model, param, effv, rtype = 'quantal', Scaled = TRUE, sav = FALSE){
	#calculate effect concentrations using associated inverse function
	if (missing(model) || missing (param)) stop('argument missing')
	if (missing(effv)) stop('error! input effv in ECx')
	if (is.vector(param)) param <- t(param)
	#effv <- sort(effv)

	if(max(effv) >= 1.0){
		rtype <- 'continuous'
		Scaled <- FALSE
	}
	
	ecx <- matrix(0, length(model), length(effv))
	
	if((rtype == 'continuous' || rtype == 'ctn') && Scaled == TRUE){
		rspnRange <- matrix(0, length(model), 2)
		rspnRange <- CEx(model, param, c(0, 1e20))
		effvAbs <- matrix(0, length(model), length(effv))
		
		for(j in seq(model)) effvAbs[j, ] <- rspnRange[j , 1] + (rspnRange[j , 2] - rspnRange[j , 1]) * effv
	}
	
	effv0 <- effv
	
	for (i in seq(model)){
		fun <- model[i]
		p <- param[i, ]
		
		if((rtype == 'continuous' || rtype == 'ctn') && Scaled == TRUE) effv0 <- effvAbs[i, ]

		ec <- switch(fun,
			'Hill' = p[1] / ((1 / effv0 - 1)^(1 / p[2])),
			'Hill_two' = p[1] * effv0 / (p[2] - effv0),
			'Hill_three' = p[1] / ((p[3] / effv0 - 1)^(1 / p[2])),			
			'Hill_four' = p[1] / (((p[3] - p[4]) / (effv0 - p[4]) - 1)^(1 / p[2])),
			'Weibull' = exp(-(-log(log(-1 / (-1 + effv0))) + p[1]) * log(10) / p[2]),
			'Weibull_three' = exp(-(-log(log(p[3] / (p[3] - effv0))) + p[1]) * log(10) / p[2]),
			'Weibull_four' = exp((log(log((-p[4] + p[3]) / (p[3] - effv0))) - p[1]) * log(10) / p[2]),
			"Logit" = exp(-log(10) * (p[1] + log(-(-1 + effv0) / (effv0))) / p[2]),
			'Logit_three' = exp(-log(10) * (p[1] + log((p[3] - effv0) / effv0)) / p[2]),
			'Logit_four' = exp(-log(10) * (p[1] + log(-(p[3] - effv0) / (p[4] - effv0))) / p[2]),
			"BCW" = exp(log(-(p[1] * p[3] - p[2] - log(-log(1 - effv0)) * p[3]) / p[2]) / p[3]),
			"BCL" = exp(log(-(p[1] * p[3] - p[2] + log(-(-1 + effv0) / effv0) * p[3]) / p[2]) / p[3]),
			"GL" = exp(-log(10) * (p[1] + log(exp(-log(effv0) / p[3]) - 1)) / p[2])
		)
		
		ecx[i, ] <- ec
	}
	
	if(rtype == 'quantal'){
		colName <- paste0('EC', effv * 100)
		
	}else if(rtype == 'continuous' || rtype == 'ctn'){
	
		if(Scaled == FALSE) colName <- paste0('EC', effv)
		
		if(Scaled == TRUE){
			colName <- paste0('EC', effv * 100)
			colNameAbs <- paste0('Abs_rspn@E', effv * 100)
			colnames(effvAbs) <- c(colNameAbs)
			
			if(is.null(rownames(param))) rownames(effvAbs) <- model else rownames(effvAbs) <- rownames(param)
		}
	}
	
	colnames(ecx) <- colName
	
	if(is.null(rownames(param))) rownames(ecx) <- model else rownames(ecx) <- rownames(param)
	
	if((rtype == 'continuous' || rtype == 'ctn') && Scaled == TRUE){
		Results <- list(ecx = ecx, effvAbs = effvAbs)
		if (sav != FALSE){
			if(sav == TRUE) {
				sav = paste("ECx_", Sys.Date(), ".txt", sep = "")
			}
		sink(sav)
		print(Results)
		sink()
	}
	}else{
		Results <- ecx
		if (sav != FALSE){
			if(sav == TRUE) {
				svfile = paste("ECx_", Sys.Date(), ".txt", sep = "")
				write.table(Results, svfile, sep = "\t", quote = F, col.names=NA)
			} else{
				write.table(Results, sav, sep = "\t", quote = F, col.names=NA)
			}
		}
	}

	return(Results)
}
