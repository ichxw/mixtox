nmECx <- function(model, param, effv, minx, gap = -1e-6, sav = FALSE){
	#calculate effect concentrations using associated inverse function
	if (missing(model) || missing (param) || missing(effv) || missing(minx)) stop('argument missing')
	if (is.vector(param)) param <- t(param)
	#effv <- sort(effv)
	
	if(min(effv) < 0){
		effv_neg <- effv[effv < 0]
		effv_pos <- effv[effv >= 0]
		len_neg <- length(effv_neg)
		len_pos <- length(effv_pos)

		ecx <- matrix(0, length(model), (len_neg * 2 + len_pos))
		
		effv_negD <- as.numeric(format(effv_neg, digit = 3))
		effv_posD <- as.numeric(format(effv_pos, digit = 3))
		left_name <- paste('ECL', effv_negD * 100, sep = '')
		right_name <- paste('ECR', effv_negD * 100, sep = '')
		if(len_pos > 0) pos_name <- paste('EC', effv_posD * 100, sep = '') else pos_name <- c()
		colName <- c(left_name, right_name, pos_name)
		
	}else{
		len_pos <- length(effv)
		len_neg <- 0
		ecx <- matrix(0, length(model), len_pos)
		effvD <- as.numeric(format(effv, digit = 3))
		colName <- paste('EC', effvD * 100, sep = '')
	}
	colnames(ecx) <- colName
	if (is.null(rownames(param)))  rownames(ecx) <- model else rownames(ecx) <- rownames(param)
	
	b <- 10000 # upper bound.
	eps <- 1e-10
	for(i in seq(model)){
		a <- minx[i]
		if(model[i] == 'Brain_Consens') f <- paste('1 - (1 + ', param[i, 1], '* xx) / (1 + exp(', param[i, 2], '*', param[i, 3],') * xx^',param[i, 2],')', sep = '')
		if(model[i] == 'BCV') f <- paste('1 -', param[i, 1], '* (1 +', param[i, 2], '* xx) / (1 + (1 + 2 *', param[i, 2], '*', param[i, 3],') * (xx /', param[i, 3],')^',param[i, 4],')', sep = '')
		if(model[i] == 'Cedergreen') f <- paste('1 - (1 +', param[i, 1], '* exp(-1 / (xx^',param[i, 2],'))) / (1 + exp(',param[i, 3], '* (log(xx) - log(', param[i, 4],')))', sep = '')
		if(model[i] == 'Beckon') f <- paste('(',param[i, 1], '+ (1 - (',param[i, 1],') / (1 + (',param[i, 2], '/ xx)^',param[i, 3],'))) / (1 + (xx /', param[i, 4],')^',param[i, 5],')', sep = '')
		if(model[i] == 'Biphasic') f <- paste(param[i, 1], '-', param[i, 1], '/ (1 + 10^((xx -', param[i, 2],') *', param[i, 3],')) + (1 -', param[i, 1],') / (1 + 10^((',param[i, 4], '- xx) *', param[i, 5],'))', sep = '')
		if(model[i] == 'Hill_five') f <- paste('1 - (1 + (', param[, 3], '- 1) / (1 + (', param[i, 1], '/xx)^', param[i, 2], ')) * (1 - 1 / (1 + (', param[i, 4], '/xx)^', param[i, 5], '))', sep = '')
		
		for(j in seq(len_pos + len_neg)){
			# the right side of minx
			fun_body <- paste(f, '-', effv[j], sep = '')
			fun = function(xx) eval(parse(text = fun_body))
			root <- uniroot(fun, c(a, b), tol = eps)$root
			ecx[i, (len_neg + j)] <- root
		}
	
		if(len_neg > 0){
			# the left side of minx
			for(k in seq(len_neg)){
				negk <- effv_neg[k]
				if(negk > gap){
					negk <- gap - 1 / 100000 # correst the response gap between response limit at extreme low concentration with theoretical limit of zero
					warning('response ', effv_neg[k], ' was set to left limit of response ', negk)
				}
				fun_body <- paste(f, '-', negk, sep = '')
				fun = function(xx) eval(parse(text = fun_body))
				root <- uniroot(fun, c(eps, a), tol = eps)$root
				ecx[i, k] <- root
			}
		}
		
	}

	if (sav != FALSE){
		if(sav == TRUE) {
			svfile = paste("nmECx_", Sys.Date(), ".txt", sep = "")
			write.table(ecx, svfile, sep = "\t", quote = F, col.names=NA)
		} else{
			write.table(ecx, sav, sep = "\t", quote = F, col.names=NA)
		}
	}

	return(ecx)
}
