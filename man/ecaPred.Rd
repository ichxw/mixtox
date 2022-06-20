\name{ecaPred}
\alias{ecaPred}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Mixture Effect Predicted by CA at Particular Concentration of a Mixture}
\description{According to the fitted concentration-response information of mixtures. The 
			concentration (e.g., ECx) that causes certain effect in the mixture will be
			calculated. ecaPred will predict how much effect will be caused at ECx according
			to concentration addition.
			The individual concentration-responses should be  
			fitted only based on the following six models: Hill, 
			Weibull, Logit, BCW, BCL, and GL.}
\usage{ecaPred(effv, sgl, mix, pctMix, sav = FALSE)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{effv}{numeric vector with single or multiple effect values (0 ~ 1).}
  \item{sgl}{A list with sgl$model and sgl$param. sgl$model is character vector of equations used
		to fit the concentration-response data of individual chemicals: Hill, Weibull, Logit, BCW,
		BCL, GL. sgl$param is numeric matrix of fitting coefficients with rownames (equations) 
		and colnames (Alpha, Beta, and Gamma). For equations with only two parameters, 
		Gamma can be set as zero or any other numeric value.
		}
  \item{mix}{A list with mix$model and mix$param. mix$model is character vector of equations used
		to fit the concentration-response data of mixtures: Hill, Weibull, Logit, BCW, BCL, GL.
		mix$param is numeric matrix of fitting coefficients of mixtures' concentration-response
		data with rownames (selected equations) and colnames (Alpha, Beta, and Gamma). For 
		equations with only two parameters, Gamma can be set as zero or any other numeric value.
		}
  \item{pctMix}{A numeric matrix, the concentration ratio (percent) of every component in the
		mixture.
		}
  \item{sav}{TRUE: save output to a default file; FALSE: output will not be saved; 
  a custom file directory: save output to the custom file directory.} 
}
\details{ecaPred calculate the effect concentrations based on the fitted concentration-response
		information of the mixture according to the input effects effv (e.g., 0.05 and 0.5). The
		concentration of individual component \eqn{c_i} is computed based on mixture's 
		ECx and the proportion of component in the mixture \eqn{p_i}. 
		Then the CA effect will be calculated based 
		on the concentration addition.
		}
\value{ A numeric matrix of predicted effects}

%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{eiaPred}}
}
\examples{
## example
# predict the CA predicted response at the concentrations that cause 5\%, 10\%, 20\%, and 50\% 
# effect of antibiotic mixtures
# each mixture contains eight components. Totally, there are 10 mixtures designed by the udcr.

sgl <- antibiotox$sgl
mix <- antibiotox$udcr.mix
pct <- antibiotox$udcr.pct
ecaPred(effv = c(0.05, 0.1, 0.20, 0.5), sgl, mix, pct)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{concentration addition}
\keyword{mixture effect}
