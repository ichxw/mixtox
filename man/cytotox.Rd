\name{cytotox}
\alias{cytotox}
\docType{data}
\title{Cytotoxicity of Heavy Metal Ions and Ionic Liquids on MCF-7}
\description{
Chemicals include four heavy metal ions: 
\eqn{{\rm{NiN}}{{\rm{O}}_{\rm{3}}}\cdot 6{{\rm{H}}_{\rm{2}}}{\rm{O}}}(Ni),
\eqn{{\rm{ZnS}}{{\rm{O}}_{\rm{4}}} \cdot {\rm{7}}{{\rm{H}}_{\rm{2}}}{\rm{O}}} (Zn), 
\eqn{{\rm{CuS}}{{\rm{O}}_{\rm{4}}} \cdot 5{{\rm{H}}_{\rm{2}}}{\rm{O}}} (Cu), and 
\eqn{{\rm{MnC}}{{\rm{l}}_{\rm{2}}} \cdot 4{{\rm{H}}_{\rm{2}}}{\rm{O}}} (Mn); 
four ionic liquids are 1-Octyl-3-methylimidazolium chloride (Omim), 
1-Dodecyl-3-methylimidazolium chloride (Dmim), 1-Ethyl-3-methylimidazolium tetrafluoroborate 
(Emim), and 1-Hexyl-3-Methylimidazolium tetrafluoroborate (Hmim). \cr
The concentration-response data include the cytotoxicity of eight compounds,
two mixtures designed by eecr, and ten mixtures designed by udcr.\cr
The fitting information of eight chemicals and a total of 12 mixtures. 
}
\usage{data(cytotox)}
\format{
  A list with concentration-response data of 8 chemicals, 2 eecr mixtures, 10 udcr mixtures,
  and associated fitting information.
  \describe{
    \item{\code{Ni$x}}{a numeric vector of test concentrations}
    \item{\code{Ni$y}}{a numeric matrix of responses}
	\item{\code{Ni$name}}{test substances}
	\item{\code{Ni$concNum}}{the number of test concentrations}
	\item{\code{Ni$tierNum}}{the number of repetitions}
	\item{\code{Ni$type}}{type of test substance: single chemicals or mixtures}
	
	\item{\code{sgl$model}}{models used to fit the concentration-response data of individual 
		chemicals}
	\item{\code{sgl$param}}{fitted coefficients of concentration-response curves 
		}
	
	\item{\code{udcr.mix$model}}{models used to fit the concentration-response data of udcr 
	mixtures}
	\item{\code{udcr.mix$param}}{fitted coefficients of the concentration-response 
	curves corresponding to udcr.mix$model
		}
	\item{\code{$udcr.pct}}{the proportion of individual chemicals in udcr mixtures }	
  }
}
\details{
Quantal responses[0, 1] are needed for curve fitting using the following six equaitons:
i.e., Weibull, Logit, Hill, BCL, GL, BCW. The following equation could transform 
continous responses to quantal ones: 
\deqn{E = \frac{{{I_0} - {I_i}}}{{{I_0}}}}
where \eqn{I_0} is the average of controls for inhibition tests or the average 
of the maximum effect for stimulation tests and \eqn{I_i} the average 
effect of the \eqn{i^{th}} treatment.
}
\source{ The cytotoxicity experiments were conducted in our lab}
\examples{
# example 1
## Retrive the toxicity data of Ni on MCF-7.
cytotox$Ni

# example 2
## Retrive the toxicity information of ten udcr mixtures on MCF-7.
cytotox$udcr.mix
}
\keyword{quantal responses}
\keyword{continuous responses}
\keyword{cytotoxicity}
