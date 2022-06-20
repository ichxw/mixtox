\name{tuneFit}
\alias{tuneFit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Find Optimal Starting values for Curve Fitting}
\description{ Curve fitting is dependent on the package minpack.lm.
This generic function first seaches optimal starting values based on 
trial and error. The the concentration response data will be fitted using the optimal 
starting values. The statistics for goodness of fit
is evaluated by the following statistics: coefficient of
determination (\eqn{R^2}), adjusted coefficient of determination (\eqn{R_{adj}^2}),
root mean squared error (RMSE), mean absolute error (MAE), Akaike information criterion (AIC),
bias-corrected Akaike information criterion(AICc), and Bayesian information criterion (BIC).
Thirteen sigmoidal models ("Hill", "Hill_two", "Hill_three", "Hill_four",
"Weibull", "Weibull_three", "Weibull_four", "Logit", "Logit_three", "Logit_four",
"BCW(Box-Cox-Weibull)", "BCL(Box-Cox-Logit)", "GL(Generalized Logit)") and four
J-shaped models ("Brain_Consens", "BCV", "Biphasic", "Hill_five")
are provided to fit the concentration-response data. 
}
\usage{
tuneFit(conc, rspn, eq = 'Weibull', effv, rtype = 'quantal', rsq = 0.6, highBar = 5000, 
		bar = 1000, sav = FALSE)}

\arguments{
  \item{conc}{a numeric vector (matrix) of experimental concentrations}
  \item{rspn}{a numeric vector (matrix) of responses corresponding to conc, it should have the
  same length (or rows and columns for matrix) as conc.}
  \item{eq}{models for curve fitting: "Hill", "Hill_two", "Hill_three", "Hill_four",
            "Weibull", "Weibull_three", "Weibull_four", "Logit", "Logit_three", "Logit_four",
            "BCW", "BCL", "GL", "Brain_Consens", "BCV", "Biphasic", "Hill_five".}
  \item{rtype}{the response type of endpoint: 'continuous' or 'quantal' data.}
  \item{effv}{numeric response to calculate effect concentration, scaled respones [0, 1].}				
  \item{rsq}{r2 below which would be ignored.}
  \item{highBar}{if the number of starting values exceeds highBar, a random sample of 
					starting values will be taken.}
  \item{bar}{the number of random samples.}
  \item{sav}{TRUE: save output to a default file; FALSE: output will not be saved; 
    a custom file directory: save output to the custom file directory.} 
}
\details{tuneFit provides high frequency trial and error approach to find appropriate starting 
	values for users. It will deploy those starting values one by one until 
	finding the right one.\cr
	Function tuneFit can also be used to fit the concentration response data for a batch of 
	chemicals and calculate corresponding effect concentration.
	}
\value{
\item{sta }{goodness of fit statistics: (\eqn{R^2}, \eqn{R_{adj}^2}, MAE, RMSE, AIC,
			AICc, and BIC)}
}
\note{
	tuneFit will load the file staval.rda which contains hundreds of starting values for each 
	of the sigmoidal and hormetic models. However, 
	those starting values are also limited. Users are encouraged to send their fitted 
	coefficients to us to extent the coverage of staval.
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
## example 1
# Fit the non-monotonic concentration-response data
# we'll get a fit with r2 of 0.740
x <- hormesis$OmimCl$x
expr <- hormesis$OmimCl$y
y <- rowMeans(expr)
tuneFit(x, y, eq = 'Biphasic')

## example 2
# Fit the non-monotonic concentration-response data
# use r2 (rsq) of 0.9, we'll get a fit with r2 of 0.989
# calcualte the effect concentration that causes 5\% inhibition
x <- hormesis$OmimCl$x
expr <- hormesis$OmimCl$y
y <- rowMeans(expr)
tuneFit(x, y, eq = 'Biphasic', effv = 0.05, rsq = 0.9)

## example 3
# Fit the concentration-response data of heavy metal Ni(2+) on MCF-7 cells.
# Calculate the concentration that causes 5\% inhibition on the growth of MCF-7

x <- cytotox$Ni$x
expr <- cytotox$Ni$y
y <- rowMeans(expr)
tuneFit(x, y, eq = 'Logit', effv = 0.05)

## example 4
# Fit the concentration-response data of Paromomycin Sulfate (PAR) on photobacteria.
# Calculate the concentrations that cause 50\% inhibition on the growth of photobacteria 

x <- antibiotox$PAR$x
expr <- antibiotox$PAR$y
y <- rowMeans(expr)
tuneFit(x, y, eq = 'Logit', effv = 0.5)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{curve fitting }
\keyword{goodness of fit}
\keyword{trial and error}
\keyword{effect concentration}
