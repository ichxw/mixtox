\name{gcaHill}
\alias{gcaHill}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Mixture Toxicity Prediction Using GCA (Hill_two)}
\description{Predict the mixture toxicity based on individual concentration-response information 
			fitted by Hill_two equation. An explicit formula for gca prediction were used 
			instead of the dichotomy algorithm in gcaPred. Three 
			optional mixture design methods are provided. 
			One is the arbitrary concentration ratio (acr) for mixture components.\cr
			Users can deign random ratios for components in the mixture. Other two 
			options are equal effect concentration ratio (eecr) and uniform design concentration 
			ratio (udcr).
			}
\usage{gcaHill(model, param, mixType, effv, refEffv = c(0.10, 0.50), rtype, sav = FALSE)}
%- maybe also 'usage' for other objects documented here.
\arguments{
     \item{model}{character vector of equation names, just Hill_two}
  \item{param}{numeric matrix of fitting coefficients with rownames (equations) and 
		colnames (Alpha, Beta).
		}
  \item{mixType}{experimental design of the mixture. acr: arbitrary concentration ratio; 
		eecr: equal effect concentration ratio; udcr: uniform design concentration ratio.}
  \item{effv}{numeric vector with single or multiple (scaled) effect values (0 ~ 1).}
  \item{refEffv}{to determine the concentration ranges for predicting effect. 
         Use scaled values (0 ~ 1).}
  \item{rtype}{ the response type of endpoint: 'continuous' or 'quantal' data.}
  \item{sav}{TRUE: save output to a default file; FALSE: output will not be saved; 
  a custom file directory: save output to the custom file directory.}

}
\details{The generalized concentration addition (GCA) model was proposed by Howard and Webster 
		to predict mixtures containing partial agonists (Howard and Webster, 2009). Empirical 
		data are used to fit concentration-response function, and then predict the mixture 
		response using the inverse function.
		\deqn{E_{mix}^{GCA} = \frac{{\sum\limits_{i = 1}^n {\frac{{{\alpha _i}{c_i}}}{{{K_i}}}}}}
		{{1 + \sum\limits_{i = 1}^n {\frac{{{c_i}}}{{{K_i}}}} }}}
		where \eqn{c_i} is the concentration of component \eqn{i} in the mixture. Parameter 
		\eqn{\alpha _i} and \eqn{K_i} are fitted coefficient of \eqn{i^{th}} component, which are 
		the same as \eqn{\beta} and \eqn{\alpha} in Hill_two equation. Right, the \eqn{\alpha _i} 
		and \eqn{K_i} are corresponding to \eqn{\beta} and \eqn{\alpha} in Hill_two equation.
		}
\value{
\item{x}{ a series of concentrations}
\item{e}{a series of effects caused by the concentrations (x) as predicted by gca}
\item{pct}{the concentration ratio (percent) of every component in the mixture}
\item{uniTab}{the uniform design table used to construct the mixture when mixType is udcr}
}
\references{
Howard, G.J., Schlezinger, J.J., Hahn, M.E., Webster, T.F., 2010. Generalized Concentration 
Addition Predicts Joint Effects of Aryl Hydrocarbon Receptor Agonists with Partial Agonists 
and Competitive Antagonists. Environ. Health Perspect. 118, 666-672.\cr
Howard, G.J., Webster, T.F., 2009. Generalized concentration addition: A method for examining 
mixtures containing partial agonists. J. Theor. Biol. 259, 469-477.\cr
Hadrup, N., Taxvig, C., Pedersen, M., Nellemann, C., Hass, U., Vinggaard, A.M., 2013. 
Concentration addition, independent action and generalized concentration 
addition models for mixture effect 
prediction of sex hormone synthesis in vitro. PLoS One 8, e70490.
}

\note{
 Only for concentration-response curves fitted by the Hill_two equation in curveFit.
}
\seealso{
\code{\link{gcaPred}}
}
\examples{
model <- c("Hill_two", "Hill_two", "Hill_two", "Hill_two")
param <- matrix(c(3.94e-5, 0.97, 0, 5.16e-4, 1.50, 0, 3.43e-6, 1.04, 0, 9.18e-6, 0.77, 0), 
				nrow = 4, ncol = 3, byrow = TRUE)
rownames(param) <- c('Ni', 'Zn', 'Cu', 'Mn')
colnames(param) <- c('Alpha', 'Beta', 'Gamma')
## example 1
# using GCA to predict the mixtures designed by equal effect concentration ratio (eecr) at
# the effect concentration of EC05 and EC50
# the eecr mixture design is based on four heavy metals (four factors).
gcaHill(model, param, mixType = "eecr", effv = c(0.05, 0.5), rtype = 'continuous')

## example 2
# using GCA to predict the mixtures designed by uniform design concentration ratio (udcr)
# the udcr mixture design is based on four heavy metals (four factors).
# Seven levels (EC05, EC10, EC15, EC20, EC25, EC30, and EC50 ) are allocated in 
# the uniform table
effv <- c(0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.50)
gcaHill(model, param, mixType = "udcr", effv, rtype = 'quantal')

## example 3
# using GCA to predict the mixtures designed by arbitrary concentration ratio (acr)
# the udcr mixture design is based on four heavy metals (four factors).
# the every component in the mixture shares exactly the same ratio (0.25) 
effv <- c(0.25, 0.25, 0.25, 0.25)
gcaHill(model, param, mixType = "acr", effv)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ generalized concentration addition }
\keyword{ Hill_two }
\keyword{ equal effect concentration ratio }
\keyword{ uniform design concentration ratio }
\keyword{ arbitrary concentration ratio }
\keyword{ uniform design table}
