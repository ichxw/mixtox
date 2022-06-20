\name{showEq}
\alias{showEq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{List Requested Equations}
\description{Show the formula of different equations upto request.
}
\usage{showEq(eq)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eq}{equation name to qury}
  }
\details{
Thirteen monotonic(sigmoidal) equations ("Hill", "Hill_two", "Hill_three", "Hill_four", 
"Weibull", "Weibull_three", "Weibull_four", "Logit", "Logit_three", "Logit_four", 
"BCW(Box-Cox-Weibull)", "BCL(Box-Cox-Logit)", "GL(Generalized Logit)") and four
non-monotonic(J-shaped) equations ("Brain_Consens", "BCV", "Biphasic", "Hill_five") 
are provided to fit concentration-response data. 

}
\value{The formula of requested equations (with abbr.) will show up.}
\references{
Scholze, M. et al. 2001. A General Best-Fit Method for Concentration-Response Curves and the
Estimation of Low-Effect Concentrations. Environmental Toxicology and Chemistry
20(2):448-457.\cr
Zhu X-W, et.al . 2013. Modeling non-monotonic dose-response relationships: Model evaluation
and hormetic quantities exploration. Ecotoxicol. Environ. Saf. 89:130-136.\cr
Howard GJ, Webster TF. 2009. Generalized concentration addition: A method for examining mixtures 
containing partial agonists. J. Theor. Biol. 259:469-477.\cr	
Spiess, A.-N., Neumeyer, N., 2010. An evaluation of R2 as an inadequate measure for nonlinear
models in pharmacological and biochemical research: A Monte Carlo approach. BMC Pharmacol.
10, 11.\cr
Di Veroli GY, Fornari C, Goldlust I, Mills G, Koh SB, Bramhall JL, et al. 2015. An automated 
fitting procedure and software for dose-response curves with multiphasic features. 
Scitific Report 5: 14701.\cr
Gryze, S. De, Langhans, I., Vandebroek, M., 2007. Using the correct intervals for prediction: A
tutorial on tolerance intervals for ordinary least-squares regression. Chemom. Intell. Lab.
Syst. 87, 147-154.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{curveFit}}
}
\examples{
# example 1
## show Weibull model
showEq('Weibull')

# example 2
## show the name of all sigmoidal models
showEq('sigmoid')
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ monotonic and non-monotonic equations }
