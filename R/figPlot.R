figPlot <- function(object, logT = TRUE, xlabel = "concentration (mol/L)", ylabel = "Response", ylimit, lgd = NULL){
	# plot the concentration-response curves
	#tiff(file = paste(root_name, "_04-12.tiff", sep = ""), res = 100)
	crcInfo <- object$crcInfo
	
	size <- dim(crcInfo)
	x <- crcInfo[, 1]
	yhat <- crcInfo[, 2]
	rspn <- crcInfo[, 3 : (size[2] - 4)]
	if(is.vector(rspn)) rspn <- as.matrix(rspn)
	PI <- crcInfo[, (size[2] - 3) : (size[2] - 2)]
	CI <- crcInfo[,  (size[2] - 1) : size[2]]
	rtype <- object$rtype
	rspnRange <- object$rspnRange
	
	if(rtype != 'continuous'){
		if(missing(ylimit)) ylimit <- c((min(rspn) * 100 -20), (max(rspn) * 100 + 20))
	}else{
		if(missing(ylimit)) ylimit <- c((rspnRange[1] - 0.1 * (rspnRange[2] - rspnRange[1])), (rspnRange[2] + 0.1 * (rspnRange[2] - rspnRange[1])))
	}
	
	par(mar=c(5,5,1,1))
	
	ifelse(logT == TRUE, xx <- log10(x), xx <- x)
	
	if(rtype != 'continuous'){
		plot(rep(xx, ncol(rspn)), rspn * 100, ylim = ylimit, xlab = xlabel, ylab = ylabel, pch = 16, cex = 2, cex.lab = 1.8, cex.axis = 1.8)
		lines(xx, yhat * 100, col = 1, lwd = 1.9)
		lines(xx, PI[, 1] * 100, col = 'blue', lwd = 1.9)
		lines(xx, PI[, 2] * 100, col = 'blue', lwd = 1.9)
		lines(xx, CI[, 1] * 100, col = 'red', lwd = 1.9)
		lines(xx, CI[, 2] * 100, col = 'red', lwd = 1.9)
	}else if(rtype == 'quantal' || rtype == 'hormesis'){
		plot(rep(xx, ncol(rspn)), rspn, ylim = ylimit, xlab = xlabel, ylab = ylabel, pch = 16, cex = 2, cex.lab = 1.8, cex.axis = 1.8)
		lines(xx, yhat, col = 1, lwd = 1.9)
		lines(xx, PI[, 1], col = 'blue', lwd = 1.9)
		lines(xx, PI[, 2], col = 'blue', lwd = 1.9)
		lines(xx, CI[, 1], col = 'red', lwd = 1.9)
		lines(xx, CI[, 2], col = 'red', lwd = 1.9)
	}
	
	if (is.null(lgd) == FALSE){
		legend('topleft', lgd, cex = 2)
	}
	#legend("topleft", inset = 0.01, root_name, box.col = 'white', cex = 1.9) 
	#dev.off()
}