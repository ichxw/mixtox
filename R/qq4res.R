qq4res <- function(object, xlabel = 'Theoretical Quantiles', ylabel = 'Residuals', lgd = NULL){
	#qqnorm(res, main = main, ylab = ylabel, cex = 2.5, pch = 16, cex.lab = 2, cex.axis = 2)
	res <- object$res
	par(mar=c(5,5,1,1))
	qqnorm(res, xlab = xlabel, ylab = ylabel, cex = 2, cex.lab = 1.8, cex.axis = 1.6)
	qqline(res, lwd = 1.8)
	
	if (is.null(lgd) == FALSE) {
		legend('topleft', lgd, cex = 2)
	}
}
