\name{gcaPred}
\alias{gcaPred}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Mixture Toxicity Prediction Using GCA (General)}
\description{Predict the mixture toxicity based on individual 
			concentration-response information. Thirteen monotonic(sigmoidal) models
			("Hill", "Hill_two", "Hill_three", "Hill_four",
		"Weibull", "Weibull_three", "Weibull_four", "Logit", "Logit_three", "Logit_four",
		"BCW(Box-Cox-Weibull)", "BCL(Box-Cox-Logit)", "GL(Generalized Logit)") are 
		incorporated to construct the GCA equation. The dichotomy 
			technique is used to solve the constructed equation. 
			Three optional mixture design methods are provided. One is the arbitrary 
			concentration ratio (acr) for mixture components. Users can arbitrarily 
			deign a random ratio for each component in the 
			mixture. Other two options are equal effect concentration 
			ratio (eecr) and uniform design concentration ratio (udcr).
			}
\usage{gcaPred(model, param, mixType, effv, refEffv, rtype, lb, ub, sav = FALSE)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{vector of equation names}
  \item{param}{numeric matrix of fitting coefficients with rownames (equation selected) and 
  colnames (ALpha, Beta, and Gamma). For equations with two parameters, Gamma can be set as 
  zero or any other numeric value.
  }
  \item{mixType}{experimental design of the mixture. acr: arbitrary concentration ratio; 
		eecr: equal effect concentration ratio; udcr: uniform design concentration ratio.
		}
  \item{effv}{numeric vector with single or multiple (scaled) effect values (0 ~ 1).}
  \item{refEffv}{to determine the concentration ranges for predicting effect. 
        Use scaled values (0 ~ 1).}
  \item{lb}{lower bound for solving constructed IA using diochotomy method (default is 1E-9).}
  \item{ub}{upper bound bound for solving constructed IA using diochotomy method(default is 9).}
  \item{rtype}{ the response type of endpoint: 'continuous' or 'quantal' data.}
  \item{sav}{TRUE: save output to a default file; FALSE: output will not be saved; 
  a custom file directory: save output to the custom file directory.}
}
\details{The generalized concentration addition (GCA) model was proposed by Howard and Webster to 
 predict mixtures containing partial agonists (Howard and Webster, 2009). 
 \deqn{\sum\limits_{i = 1}^n {\frac{{{c_i}}}{{f_i^{ - 1}\left( E \right)}}}  = 1}
 Empirical data were used to fit concentration-response function, and then predict the mixture 
 response using the inverse function. Previous studies used Hill_two function to fit individual 
 concentration response curves in the GCA prediction (Hadrup et al., 2013; Howard et al., 2010). 
 Here, we incorporated thirteen functions to construct
 the GCA equation and the dichotomy technique is used to solve the constructed equation. 
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
mixtures containing partial agonists. J. Theor. Biol. 259, 469-477.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
\code{\link{gcaHill}}
}
\examples{
## example 1
# using GCA to predict the mixtures designed by equal effect concentration ratio (eecr) at the 
# effect concentration of EC50
# the eecr mixture design is based on seven antibiotics(seven factors).
model <- antibiotox$sgl$model
param <- antibiotox$sgl$param
refEffv <- c(0.1, 0.50, 0.80)
gcaPred(model, param, mixType = "eecr", effv = 0.5, refEffv, rtype = 'quantal')

## example 2
# using GCA to predict the mixtures designed by uniform design concentration ratio (udcr)
# the udcr mixture design is based on 2 antibiotics(2 factors) and 
# three levels (EC05, EC20, and EC50 )
model <- antibiotox$sgl$model[1 : 2]
param <- antibiotox$sgl$param[1 : 2, ]
effv <- c(0.05, 0.20, 0.50)
refEffv <- c(0.1, 0.80)
gcaPred(model, param, mixType = "udcr", effv, refEffv, rtype = 'quantal')

## example 3
# using GCA to predict the mixtures designed by arbitrary concentration ratio (acr)
# the udcr mixture design is based on 2 heavy metals (2 factors).
# the every component in the mixture shares exactly the same ratio (0.5) 
model <- cytotox$sgl$model[1 : 2]
param <- cytotox$sgl$param[1 : 2, ]
effv <- c(0.5, 0.5)
refEffv <- c(0.1, 0.80)
gcaPred(model, param, mixType = "acr", effv, refEffv, rtype = 'quantal')
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ generalized concentration addition }
\keyword{ equal effect concentration ratio }
\keyword{ uniform design concentration ratio }
\keyword{ arbitrary concentration ratio }
\keyword{ uniform design table}
