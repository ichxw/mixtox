NOEC <- function(x, rspn, blankC = FALSE, sigLev = 0.05, alternertive = 'B', sav = FALSE){
	# NOEC and LOEC calculation using Dunnett's test
	## Dunnett, C.W., 1964. New tables for multiple comparisons with controls. Biometrics 30, 482-491
	## Q: One dataset has four blank controls (C1, C2, C3, C4) and one treatment has three replicates (T1, T2, T3), 
	##   another treatment has five replicates (R1, R2, R3, R4, R5), how to arrange the response matrix (rspn)?
	## A: Label the missing values as NA, the response matrix (rspn) can be arranged as follows:
	##
	##		C1	C2	C3	C4	NA
	##		T1	T2	T3 NA	NA
	##   	R1	R2	R3	R4	R5
	## The adjustation of critical value for the unequal variances or unequal number of control and replicates is skipped in this program.
	## We expect the nubmer of controls and replicates is equal. 
	## non-observed effect concentration
	## least observed effect concentration
	##
	
	getDTcv <- function(n, DF, alternertive, sigLev){
		#load('DTcv.rda')
		d.f. <- c(5 : 30, 40, 50, 60, 80, 100, 120, 200, 1000, 3000)
		#rowLength <- length(d.f.)
		K <- 30
		D <- length(which(DF >= d.f.))
		rowPosition <- (D - 1) * K + n
		
		if(alternertive == 'B') {
			if(sigLev == 0.01) colPosition <- 3
			if(sigLev == 0.05) colPosition <- 4
			if(sigLev == 0.10) colPosition <- 5
		}else {
			if(sigLev == 0.01) colPosition <- 6
			if(sigLev == 0.05) colPosition <- 7
			if(sigLev == 0.10) colPosition <- 8		
		}
		
		criticalValue <- DTcv[rowPosition, colPosition]
		return(criticalValue)
	
	}
	
	if(is.vector(rspn)){
		stop("y should a matrix")
		
	}else if(is.matrix(rspn)){
		size <- dim(rspn)
		nrep <- size[2]
		if(nrep < 3) stop("y should be a response matrix with at least three replicates")
	}
	
	if(is.vector(x)) {
		if(sum(diff(x) < 0) > 0) stop('error! concentrations or levels in an ascending order')
	}else {
		x <- as.vector(x)
	}
	if(blankC != TRUE && blankC != FALSE) blankC = FALSE
	if(alternertive != 'B' && alternertive != 'U')  alternertive = 'B'
	if(sigLev != 0.01 && sigLev != 0.05 && sigLev != 0.10) sigLev == 0.05
	
	n <- nrow(rspn)
	m <- ncol(rspn)
	
	if(blankC == TRUE) {
		# blank controls are in the first row of the response matrix(rspn)
		# The first concentration is x is the blank control treatment
		blankControl <- rspn[1, ]
		x <- x[-1]
		n <- n - 1
	}else {
		blankControl <- rep(0, m)
		rspn <- rbind(blankControl, rspn)
	}
	
	SST <- 0
	SSB <- 0
	Tj <- rep(0, n + 1)
	nj <- rep(0, n + 1)
	sumSq <- 0  # the sum of all (treatment effect)^2, sum(sum(xji^2))
	C <- 0
	DT <- rep(0, n)
	for(j in seq_len(n + 1)) {
		treat <- rspn[j, !is.na(rspn[j,])]
		Tj[j] <- sum(treat)
		nj[j] <- length(treat)
		sumSq <- sumSq + sum(treat^2)
	}
	C <- sum(Tj)
	N <- sum(nj)
	DF <- N -  (n + 1)
	
	SSB <- sum(Tj^2 / nj) - C^2 / N
	SST <- sumSq - C^2 / N
	SSW <- SST - SSB # variance of treatment
	SW <- sqrt(SSW / DF) # standard error of treatment
	
	for(j in seq_len(n)) {
		DT[j] <- (Tj[j + 1] / nj[j + 1] - sum(blankControl, na.rm = TRUE) / nj[1]) / (SW * sqrt((1 / nj[j + 1]) + (1 / nj[1])))
	}
	
	#DTcv <- qf(1 - sigLev / 2, (n + 1) * (m - 1), n)
	
	DTcv <- getDTcv(n, DF, alternertive, sigLev)
	noecSign <- sign(abs(DT) - DTcv)

	idx <- which(noecSign == 1)
	if(length(idx) == 0) { 
		noec <- x[length(x)]
		loec = NULL
	}else if (length(idx) == n) {
		noec <- NULL
		loec <- x[1]
	}else {
		noec <- x[idx[1] - 1]
		loec <- x[idx[1]]
	}

	mat <- cbind(x, DT, DTcv, noecSign)
	colnames(mat) <- c('C/Level', 't', 'critical_value', 'sign')
	Results <- list(mat = mat, noec = noec, loec = loec, sigLev = sigLev, DF = c(n, DF))

	if (sav != FALSE){
		if(sav == TRUE) {
			sav = paste("NOEC_", Sys.Date(), ".txt", sep = "")
		}
		sink(sav)
		print(Results)
		sink()
	}

	return(Results)
}
