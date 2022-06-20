\name{curveFit}
\alias{curveFit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Curve Fitting}
\description{
Thirteen monotonic(sigmoidal) models ("Hill", "Hill_two", "Hill_three", "Hill_four",
"Weibull", "Weibull_three", "Weibull_four", "Logit", "Logit_three", "Logit_four",
"BCW(Box-Cox-Weibull)", "BCL(Box-Cox-Logit)", "GL(Generalized Logit)") and four
non-monotonic(J-shaped) models ("Brain_Consens", "BCV", "Biphasic", "Hill_five")
are provided to fit dose-response data. The goodness of fit of a model
is evaluated by the following statistics: coefficient of
determination (\eqn{R^2}), adjusted coefficient of determination (\eqn{R_{adj}^2}),
root mean squared error (RMSE), mean absolute error (MAE), Akaike information criterion (AIC),
bias-corrected Akaike information criterion(AICc), and Bayesian information criterion (BIC).
}
\usage{
curveFit(x, rspn, eq , param, effv, rtype = 'quantal', sigLev = 0.05, sav = FALSE, ...)}

\arguments{
\item{x}{a numeric vector of experimental concentration.}
\item{rspn}{a numeric matrix of experimental responses with one or more replicates.}
\item{eq}{equation used for curve fitting: "Hill", "Hill_two", "Hill_three", "Hill_four",
           "Weibull", "Weibull_three", "Weibull_four", "Logit", "Logit_three", "Logit_four",
           "BCW", "BCL", "GL", "Brain_Consens", "BCV", "Biphasic", "Hill_five".}
\item{param}{a vector of starting parameters. Use tuneFit to get the starting values.}
\item{effv}{a numeric vector of responses for the calculation of effect concentrations.
			Minus values(e.g., -5\%) are permited only in the condition of 'hormesis' 
			dose-responses. 
			Relative values(e.g., 5\%, 10\%) in the condition of 'continuous' dose-responses.
			}
\item{rtype}{three dose-response types: 'quantal', 'continuous', 'hormesis'. 
			Default is 'quantal'. 
			'quantal': dose-responses with lower limit fixed at 0 and higher limit at 1 (100\%). 
			'continuous': dose-responses with no fixed lower or higher limits.
			'hormesis': non-monotonic J or U-shaped dose-responses with lower limit fixed at 0 
			and higher limit at 1 (100\%).
				}
\item{sigLev}{the significant level for confidence intervals and Dunnett\'s test. Default
        is 0.05.}
\item{sav}{TRUE: save output to a default file; FALSE: output will not be saved; 
  a custom file directory: save output to the custom file directory.} 
\item{...}{other arguments passed to nlsLM in minpack.lm.}	
}
\details{Curve fitting is dependent on the package minpack.lm
        (http://cran.r-project.org/web/packages/minpack.lm/index.html).\cr
		Monotonic(sigmoidal) equations are as follows:\cr
		Hill: \deqn{E = 1/\left( {1 + {{\left( {\alpha /c} \right)}^\beta }} \right)}
		Hill_two: \deqn{E = \beta c/\left( {\alpha  + c} \right)}
		Hill_three: \deqn{E = \gamma /\left( {1 + {{\left( {\alpha /c} \right)}^\beta }} \right)}
		Hill_four: \deqn{E = \delta  + \left( {\gamma  - \delta } \right)/\left( {1 +
							{{\left( {\alpha /c} \right)}^\beta }} \right)}
		where \eqn{\alpha} = EC50, \eqn{\beta} = H (Hill coefficient), \eqn{\gamma} = Top,
		and \eqn{\delta} = Bottom \cr

		Weibull: \deqn{E = 1 - \exp ( - \exp (\alpha  + \beta \log (c)))}
		Weibull_three: 	\deqn{E = \gamma \left( {1 - \exp \left( { - \exp \left( {\alpha  +
		                       \beta \log \left( c \right)} \right)} \right)} \right)}
		Weibull_four: \deqn{E = \gamma  + \left( {\delta  - \gamma } \right)\exp \left(
		                   { - \exp \left({\alpha  + \beta \log \left( c \right)} \right)} \right)}
		Logit: \deqn{E = {(1 + \exp ( - \alpha  - \beta \log (c)))^{ - 1}}}
		Logit_three: \deqn{E = \gamma /\left( {1 + \exp \left( {\left( { - \alpha } \right)
		                - \beta \log \left( c \right)} \right)} \right)}
		Logit_four: \deqn{E = \delta  + \left( {\gamma  - \delta } \right)/\left( {1 + \exp \left
		            ( {\left( { - \alpha } \right) - \beta \log \left( c \right)} \right)} \right)}

		where \eqn{\alpha} is the location parameter and \eqn{\beta} slope parameter.
		\eqn{\gamma} = Top, and \eqn{\delta} = Bottom \cr

		BCW: \deqn{E = 1 - \exp \left( { - \exp \left( {\alpha  + \beta \left(
				{\frac{{{c^\gamma } - 1}}{\gamma }} \right)} \right)} \right)}
		BCL: \deqn{E = {(1 + \exp ( - \alpha  - \beta (({c^\gamma } - 1)/\gamma )))^{ - 1}}}
		GL: \deqn{E = 1/{(1 + \exp ( - \alpha  - \beta \log (c)))^\gamma }}

		Non-monotonic(J-shaped) models:\cr
		Hill_five: \deqn{E = 1 - \left( {1 + \left( {\gamma  - 1} \right)/\left( {1 + {{\left( {\alpha /c}
					\right)}^\beta }} \right)} \right)\left( {1 - 1/\left( {1 + {{\left( {\delta /c} \right)}
					^\varepsilon }} \right)} \right)}

		Brain_Consens: \deqn{E = 1 - \left( {1 + \alpha c} \right)/\left( {1 + \exp \left(
		                      {\beta \gamma } \right){c^\beta }} \right)}
		where \eqn{\alpha} is the initial rate of increase at low concentration, \eqn{\beta} 
		the way in which 
		response decreases with concentration, and \eqn{\gamma} no simple interpretation.\cr

		BCV: \deqn{E = 1 - \alpha \left( {1 + \beta c} \right)/\left( {1 + \left( {1 + 2\beta
					\gamma } \right){{\left( {c/\gamma } \right)}^\delta }} \right)}
		where \eqn{\alpha} is untreated control, \eqn{\beta} the initial rate of 
		increase at low concentration, \eqn{\gamma} the concentration cause 50\% inhibition, 
		and \eqn{\delta} no simple interpretation.\cr

		Cedergreen: \deqn{E = 1 - \left( {1 + \alpha \exp \left( { - 1/\left( {{c^\beta }}
						\right)} \right)} \right)/\left( {1 + \exp \left( {\gamma \left({\ln
						\left( c \right) - \ln \left( \delta  \right)} \right)} \right)} \right)}
		where \eqn{\alpha} the initial rate of increase at low concentration, \eqn{\beta} 
		the rate of	hormetic effect manifests itself, 
		\eqn{\gamma} the steepness of the curve after 
		maximum hormetic effect, and \eqn{\delta} the lower bound on the EC50 level.\cr

		Beckon: \deqn{E = \left( {\alpha  + \left( {1 - \alpha } \right)/\left( {1 + {{\left(
						{\beta /c} \right)}^\gamma }} \right)} \right)/\left( {1 + {{\left(
						{c/\delta } \right)}^\varepsilon }} \right)}
		where \eqn{\alpha} is the minimum effect that would be approached by 
		the downslope in the absence 
		of the upslope, \eqn{\beta} the concentration at the midpoint of the falling slope, 
		\eqn{\gamma} the steepness of the rising(positive) slope, \eqn{\delta} the concentration 
		at the midpoint of the rising slope, and \eqn{\epsilon} the steepness of the 
		falling(negative) slope.\cr

		Biphasic: \deqn{E = \alpha  - \alpha /\left( {1 + {{10}^{\left( {\left( {c - \beta }
						\right)\gamma } \right)}}} \right) + \left( {1 - \alpha } \right)/\left
						( {1 + {{10}^{\left( {\left( {\delta  - c} \right)\varepsilon } \right)}}}
						\right)}
		where \eqn{\alpha} is the minimum effect that would be approached by the 
		downslope in the absence 
		of the upslope, \eqn{\beta} the concentration at the midpoint of the falling slope, 
		\eqn{\gamma} the steepness of the rising(positive) slope, \eqn{\delta} the concentration 
		at the midpoint of the
		rising slope, and \eqn{\epsilon} the steepness of the falling(negative) slope.\cr
		In all, \eqn{E} represents effect and \eqn{c} represents concentration.
		
		
	}
\value{
\item{fitInfo }{curve fitting information.}
\item{eq }{equation used in curve fitting.}
\item{p }{fitted parameters.}
\item{res }{residual.}
\item{sta }{goodness of fit.}
\item{crcInfo }{a numeric matrix with the experimental concentration (x), predicted and
		experimental responses,
		experimental responses, lower and upper bounds of (non-simultaneous) prediction intervals
		(PI.low and PI.up), and lower and upper bounds of (non-simultaneous) confidence intervals
		(CI.low and CI.up).
		}
\item{ecx}{effect concentrations only if effv is provided.}
\item{effvAbs}{Absolute effects corresponding to effv only in the condition of 'continuous' 
		dose-responses.}
\item{rtype }{dose-response type.}
\item{rspnRange}{response range. The lower limit is the response at extremely low dose.
				The higher limit is the response at infinite high dose.
				}
\item{minx }{concentration to induce the maximum stimulation for 'continuous' dose-response}
\item{miny }{the maximum stimulation for 'continuous' data.}
}
\references{
Scholze, M. et al. 2001. A General Best-Fit Method for dose-response Curves and the
Estimation of Low-Effect Concentrations. Environmental Toxicology and Chemistry
20(2):448-457.\cr
Zhu X-W, et.al. 2013. Modeling non-monotonic dose-response relationships: Model evaluation
 and hormetic quantities exploration. Ecotoxicol. Environ. Saf. 89:130-136.\cr
Howard GJ, Webster TF. 2009. Generalized concentration addition: A method for examining mixtures 
containing partial agonists. J. Theor. Biol. 259:469-477.\cr	
Spiess, A.-N., Neumeyer, N., 2010. An evaluation of R2 as an inadequate measure for nonlinear
models in pharmacological and biochemical research: A Monte Carlo approach. BMC Pharmacol.
10, 11.\cr
Huet, S., Bouvier, A., Poursat, M.-A., Jolivet, E., 2004. Statistical tools for nonlinear 
regression: a practical guide with S-PLUS and R examples. Springer Science & Business Media.\cr
Gryze, S. De, Langhans, I., Vandebroek, M., 2007. Using the correct intervals for prediction: A
tutorial on tolerance intervals for ordinary least-squares regression. Chemom. Intell. Lab.
Syst. 87, 147-154.
}
\note{
	tuneFit is recommended to find the starting values.
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
## example 1
# Fit hormesis dose-response data.
# Calculate the concentrations that cause 5\% of 50\% inhibition.
x <- hormesis$OmimCl$x
rspn <- hormesis$OmimCl$y
curveFit(x, rspn, eq = 'Biphasic', param = c(-0.34, 0.001, 884, 0.01, 128), 
			effv = 0.5, rtype = 'hormesis')

x <- hormesis$HmimCl$x
rspn <- hormesis$HmimCl$y
curveFit(x, rspn, eq = 'Biphasic', param = c(-0.59, 0.001, 160,0.05, 19),  
			effv = c(0.05, 0.5), rtype = 'hormesis')

x <- hormesis$ACN$x
rspn <- hormesis$ACN$y
curveFit(x, rspn, eq = 'Brain_Consens', param = c(2.5, 2.8, 0.6, 2.44),  
			effv = c(0.05, 0.5), rtype = 'hormesis')

x <- hormesis$Acetone$x
rspn <- hormesis$Acetone$y
curveFit(x, rspn, eq = 'BCV', param = c(1.0, 3.8, 0.6, 2.44),  effv = c(0.05, 0.5), 
			rtype = 'hormesis')

## example 2
# Fit quantal dose-responses: the inhibition of heavy metal Ni(2+) on the growth of MCF-7 cells.
# Calculate the concentrations that cause 5\% and 50\% inhibition. 
x <- cytotox$Ni$x
rspn <- cytotox$Ni$y
curveFit(x, rspn, eq = 'Logit', param = c(12, 3), effv = c(0.05, 0.5), rtype = 'quantal')

## example 3
# Fit quantal dose-responses: the inhibition effect of Paromomycin Sulfate (PAR) on photobacteria.
# Calculate the concentrations that cause 5\% and 50\% inhibition.
x <- antibiotox$PAR$x
rspn <- antibiotox$PAR$y
curveFit(x, rspn, eq = 'Logit', param = c(26, 4), effv = c(0.05, 0.5))
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{curve fitting }
\keyword{dose-response curve}
\keyword{goodness of fit}
\keyword{quantal dose-reponses}
\keyword{continuous dose-reponses}
\keyword{hormesis dose-reponses}
\keyword{response range}
\keyword{effect concentration}
