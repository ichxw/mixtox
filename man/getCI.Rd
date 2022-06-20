\name{getCI}
\alias{getCI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculating Confidence Intervals }
\description{calculating non-simultanous confidence intervals and prediction intervals}
\usage{getCI(object, effv, Scaled = TRUE, sigLev = 0.05, sav = FALSE)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ object of class curveFit.}
  \item{effv}{numeric matrix of experimental responses with at least three replicates.}
  \item{Scaled}{indicating if effv was scaled or not(TRUE/FALSE) in continuous dose-response 
                (rtype = 'continuous')}
  \item{sigLev}{significance level(default is 0.05).}
  \item{sav}{TRUE: save output to a default file; FALSE: output will not be saved; 
    a custom file directory: save output to the custom file directory.} 
}
\details{ 
The Delta method (Dybowski et al, 2001) is used to construct confidence
intervals for predicted responses. 
}
\value{
\item{xmat }{effect concentration(s) and corresponding CIs and PIs}
\item{emat }{effect(s) and and corresponding CIs and PIs}
\item{sav}{TRUE: save output to a default file; FALSE: output will not be saved; 
  a custom file directory: save output to the custom file directory.} 
}
\references{ Zhu, X.-W. and Chen, J.-Y. (2016). mixtox: An R Package for Mixture 
            Toxicity Assessment. R Journal, 8(2).\cr
Dybowski, R. and Gant, V. (2001). Clinical applications of artificial neural networks. 
Cambridge University Press, Cambridge.\cr
Gryze, S. De, Langhans, I., and Vandebroek, M. (2007). Using the correct intervals for 
prediction: A tutorial on tolerance intervals for ordinary least-squares regression. 
Chemom. Intell. Lab. Syst., 87, 147-154.
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\examples{
## example 1
x <- cytotox$Ni$x
rspn <- cytotox$Ni$y
obj <- curveFit(x, rspn, eq = 'Logit', param = c(12, 3), effv = c(0.05, 0.5), rtype = 'quantal')
getCI(obj, effv = c(0.05, 0.50))
}
\keyword{non-simultanous confidence intervals}
\keyword{non-simultanous prediction intervals}
