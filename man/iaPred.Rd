\name{iaPred}
\alias{iaPred}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Mixture Toxicity Prediction Based on Independent Action}
\description{Predict the mixture toxicity based on individual concentration-response 
           information fitted only based on the following six models: Hill, 
		   Weibull, Logit, BCW, BCL, and GL. Three optional mixture design methods 
			are provided:(1)arbitrary concentration ratio (acr), users can set 
			an arbitary proportion for each component in a mixture; (2)equal effect 
			concentration ratio (eecr); (3) uniform design concentration ratio (udcr).
			}
\usage{iaPred(model, param, mixType, effv, effPoints, lb = 1e-9, ub = 6, sav = FALSE)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{ character vector of models: Hill, Weibull, Logit, BCW, BCL, GL}
  \item{param}{numeric matrix of fitting coefficients with row names (selected equations) and 
		column names (Alpha, Beta, and Gamma). For equations with two parameters, Gamma can be 
		set to zero or any other numeric value.
		}
  \item{mixType}{experimental design of the mixture. 
		acr: arbitrary concentration ratio;  
		eecr: equal effect concentration ratio;
		udcr: uniform design concentration ratio.
		}
  \item{effv}{numeric vector with single or multiple effect values}
  \item{effPoints}{numeric vector [0 ~ 1] to predict effect concentrations.}
  \item{lb}{lower bound for solving constructed IA using diochotomy method.}
  \item{ub}{upper bound bound for solving constructed IA using diochotomy method.}
  \item{sav}{TRUE: save output to a default file; FALSE: output will not be saved; 
  a custom file directory: save output to the custom file directory.} 
}

\details{IA is designed for mixtures of chemicals that have distinct 
		mechanisms of action. The IA model is commonly defined as:
		\deqn{{\rm{E}}({c_{mix}}) = 1 - (1 - {\rm{E}}({c_1}))(1 - {\rm{E}}({c_2})) \cdot  
		\cdot  \cdot (1 - {\rm{E}}({c_n})) = 1 - \prod\limits_{i = 1}^n {\left( {1 - 
		{\rm{E}}({c_i})} \right)}}
		where \eqn{E\left( {{c_{mix}}} \right)} is the overall effect caused by \eqn{c_{mix}}, 
		and \eqn{E({c_i})} is the effect elicited by \eqn{c_i} when applied individually. 
		For a fitted function \eqn{f_i} based on the concentration-response data of the 
		\eqn{i^{th}} component, \eqn{E({c_i})} is equal to \eqn{{f_i}({c_i})}.
		When \eqn{E({c_{mix}}) = x}, the equation can be expressed as: 
		\deqn{x\%  = 1 - \prod\limits_{i = 1}^n {\left( {1 - {f_i}\left( {{p_i}\left( 
		{E{C_{x,mix}}} \right)} \right)} \right)}}
		This equation can be used to predict the combined effects of mixture-components based
		on IA. The dichotomy technique is used to solve the constructed equation.
		}
\value{
\item{ia}{ a series of effect concentrations predicted by IA}
\item{e}{a series of effects (effPoints) associated with the effect concentrations in ia}
\item{pct}{the proportion of every component in a mixture}
\item{uniTab}{the uniform design table to construct the mixture when mixType is 'udcr'}
}
\references{
Liang, Yi-zeng, Kai-tai Fang, and Qing-song Xu. 2001. Uniform Design and Its Applications in 
Chemistry and Chemical Engineering. Chemometrics and Intelligent Laboratory Systems 
58(1):43-57.\cr
Backhaus, T., Faust, M., 2012. Predictive environmental risk assessment of chemical mixtures: A 
conceptual framework. Environmental Science and Technology. 46, 2564-2573.
}

\note{
Note that effv is dependent on the mixType.\cr
If the mixType is acr, effv is supposed to be a vector with the ratio of each components and its length 
should be the same as the number of components. For instance, effv can be c(4, 1) if one wants to mix two 
compounds with a ratio of 4:1. One can also set effv as c(80, 20) or c(0.8, 0.2) or any other values 
with the same ratio as 4:1. The program will convert the ratio to a range of 0 to 100\%.\cr
If the mixType is eecr, effv is supposed to be the effect values. For instance, effv = c(0.5) means one mixture 
will be prepared and each of its components will be mixed in a ratio equal to EC50. More than one eecr mixtures 
can be achieved by assigning more effect values to effv, e.g., effv = c(0.05, 0.10, 0.5) \cr
If the mixType is udcr, elements in effv are levels in a uniform table, the length of effv is the same as 
the number of runs. the number of runs should be consistent with the number of components in the uniform design.
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{caPred}}
}
\examples{
# data(cytotox)

## example 1
# using IA to predict the mixtures designed by equal effect concentration ratio (eecr) at the 
# effect concentration of EC05 and EC50
# the eecr mixture design is based on four heavy metals and four ion liquids(eight factors).
model <- cytotox$sgl$model
param <- cytotox$sgl$param
iaPred(model, param, mixType = "eecr", effv = c(0.05, 0.5))

## example 2
# using IA to predict the mixtures designed by uniform design concentration ratio (udcr)
# the udcr mixture design is based on seven antibiotics (seven factors).
# five levels (EC05, EC10, EC20, EC30, and EC50 ) are allocated in the uniform table using the 
# pseudo-level technique (Liang et al., 2001)
model <- antibiotox$sgl$model
param <- antibiotox$sgl$param
effv <- c(0.05, 0.05, 0.10, 0.10, 0.20, 0.20, 0.30, 0.30, 0.50, 0.50)
iaPred(model, param, mixType = "udcr", effv)

## example 3
# using IA to predict the mixtures designed by arbitrary concentration ratio (acr)
# the udcr mixture design is based on four antibiotics (four factors).
# the every component in the mixture shares exactly the same ratio (0.25) 
model <- antibiotox$sgl$model[1 : 4]
param <- antibiotox$sgl$param[1 : 4, ]
effv <- c(0.25, 0.25, 0.25, 0.25)
iaPred(model, param, mixType = "acr", effv)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ independent action }
\keyword{ equal effect concentration ratio }
\keyword{ uniform design concentration ratio }
\keyword{ arbitrary concentration ratio }
\keyword{ uniform design table}
