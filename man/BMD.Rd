\name{BMD}
\alias{BMD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculating benchmark dose (BMD) and lower limit of benchmark dose (BMDL)}
\description{Calculation of BMD and BMDL for both quantal and continuous dose responses. Six 
  2- or 3-paramter models ('Hill', 'Weibull', 'Logit', 'Weibull_three', 'Hill_three', 
  'Logit_three') were employed for quantal dose responses. Three 4-parameter 
  models ('Weibull_four', 'Logit_four', 'Hill_four') were employed for continuous data. 
}
\usage{BMD(object, bmr = 0.10, backg = 0, def = 'additional', eq = 'as.is', 
         sigLev = 0.05, ci = 'CI', sav = FALSE)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ object of class curveFit.}
  \item{bmr}{ numeric vector of bench mark response levels for which to calculate benchmark
             doses (should be between 0 and 1)}
 \item{backg}{numeric value specifying the background level (defaults to 0)}
 \item{def}{character string specifiying the definition of the benchmark dose to use in the
           calculations. "excess" and "additional" are for binomial response whereas "relative" 
           and "hybrid" (additive hybrid) are for continuous response.}
 \item{eq}{default list of equations.}
 \item{sigLev}{the significance level for Dunnett's test. The default is 0.05.}
 \item{ci}{types of confidence intervals (CI or PI).CI: non-simultaneous confidence intervals;
			PI: non-simultaneous prediction intervals.}  
  \item{sav}{TRUE: save output to a default file; FALSE: output will not be saved; 
  a custom file directory: save output to the custom file directory.} 
}
\details{ Selecting the Benchmark Response Level (BMR) (https://www.epa.gov/bmds/
          benchmark-dose-bmd-methods#bmr).BMR is usually set as 0.10.
}
\value{
\item{bmds }{values of BMDL, BMD, and BMDU.}
}
\references{
Benchmark Dose Technical Guidance, Risk Assessment Forum, U.S. Environmental Protection Agency, 
Washington, DC 20460, EPA/100/R-12/001, June 2012\cr
Kennyp (2002). Critical Issues in Benchmark Calculations from Continuous Data. Crit. 
Rev. Toxicol., 32, 133-153.
}
\note{ three default equations (Hill, Weibull, and Logit) were used to calculate BMD for
 quantal dose response. Three default equations('Weibull_four', 'Logit_four', 'Hill_four') 
 were used to calculate BMD for continuous dose response. BMD calculation is only availabale for
 monotonic dose reponses in this version.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
## example 1
# calcualte the BMD of heavy metal Ni(2+) on the MCF-7 cells 
x <- cytotox$Ni$x
rspn <- cytotox$Ni$y
obj <- curveFit(x, rspn, eq = 'Logit', param = c(12, 3), effv = c(0.05, 0.5), rtype = 'quantal')
BMD(obj, bmr = 0.10, backg = 0, def = 'additional', eq = 'default', sigLev = 0.05, ci = 'CI')
}
\keyword{BMD}
\keyword{BMDL}
