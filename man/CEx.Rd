\name{CEx}
\alias{CEx}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Effect Calculation for All Ninteen Curves}
\description{Calculating responses at particular concentrations.
}
\usage{CEx(model, param, conc, sav = FALSE)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{a character vector of equation names}
  \item{param}{a numeric matrix of fitting coefficients with rownames (equation selected) and 
  colnames (ALpha, Beta, and Gamma). For equations with two parameters, Gamma can be set as 
  zero or any other numeric value.
  }
  \item{conc}{a numeric vector with single or multiple concentrations.}
  \item{sav}{TRUE: save output to a default file; FALSE: output will not be saved; 
  a custom file directory: save output to the custom file directory.} 
}
\details{Responses will be calculated with provided equations (model), associated 
fitting parameters (param), and concentrations.
}
\value{ 
\item{effv}{a numeric vector of effect(s)}
}
\references{
Zhu X-W, et.al. 2013. Modeling non-monotonic dose-response relationships: Model evaluation 
and hormetic quantities exploration. Ecotoxicol. Environ. Saf. 89:130-136.\cr
Hill equation (biochemistry) http://en.wikipedia.org/wiki/Hill_equation_(biochemistry)\cr
Scholze, M. et al. 2001. A General Best-Fit Method for Concentration-Response Curves and the
 Estimation of Low-Effect Concentrations. Environmental Toxicology and Chemistry 20(2):448-457.
 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
## example 1
# calculate the responses of hormesis curves at the concentration of 0.1 and 0.02 mol/L
model <- hormesis$sgl$model
param <- hormesis$sgl$param
CEx(model, param, conc = c(0.1, 0.02)) 

## example 2
# calculate the effect caused by four heavy metals and four ionic liquids at the concentration of
# 0.00001 and 0.00002 mol/L on the MCF-7 cells
model <- cytotox$sgl$model
param <- cytotox$sgl$param
CEx(model, param, conc = c(0.00001, 0.00002)) 

## example 3
# calculate the response ranges
model <- hormesis$sgl$model
param <- hormesis$sgl$param
CEx(model, param, conc = c(0, 1e20))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ effect concentration }
\keyword{inverse function}
