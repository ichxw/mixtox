\name{antibiotox}
\alias{antibiotox}
\docType{data}
\title{Toxicity of Seven Antibiotics on Photobacteria}
\description{
Seven antibiotics are: Neomycin sulfate(NEO), streptomycin sulfate(STR), 
kanamycin sulfate(KAN), spectinomycin dihydrochloridehydrate(SPE), paromomycin 
sulfate(PAR), dihydrostreptomycin sesquisulfate hydrate(DIH), and gentamycin 
sulfate(GEN). Their toxicity on photobacteria Q67 were tested using 
microplate toxicity analysis.\cr
The concentration-responses include the toxicity of seven antibiotics,
 two eecr mixtures, and ten mixtures designed by udcr. The curve fitting 
 information of seven antibiotics and a total of 12 mixtures. 
}
\usage{data(antibiotox)}
\format{
  A list with concentration-response data of 7 antibiotics, 12 mixtures, and associated
  fitting information.
  \describe{
    \item{\code{PAR$x}}{a numeric vector of concentrations}
    \item{\code{PAR$y}}{a numeric matrix of responses}
	\item{\code{PAR$name}}{name of test substance}
	\item{\code{PAR$concNum}}{the number of concentrations}
	\item{\code{PAR$tierNum}}{the number of repetitions}
	\item{\code{PAR$type}}{type of test substance}
	
	\item{\code{sgl$model}}{models used to fit the concentration-response data
		}
	\item{\code{sgl$param}}{fitted coefficients of concentration-response curves
		}
	
	\item{\code{udcr.mix$model}}{models used to fit the concentration-response data of udcr mixtures on 
		photobacteria
		}
	\item{\code{udcr.mix$param}}{fitted coefficients of concentration-response curves
		}
	\item{\code{udcr.pct}}{the percentage of individual chemicals in the udcr mixtures}
	
	\item{\code{eecr.mix$model}}{models used to fit the concentration-response data of eecr mixtures on 
		photobacteria
		}
	\item{\code{eecr.mix$param}}{fitted coefficients of concentration-response curves
		}
	\item{\code{eecr.pct}}{the proportion of individual chemicals in the eecr mixtures}
	
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
\examples{
# example 1
## Retrive the toxicity information of PAR on photobacteria.
antibiotox$PAR
# example 2
## Retrive the toxicity information of two eecr mixtures on photobacteria.
antibiotox$eecr.mix
}
\keyword{quantal responses}
\keyword{continuous responses}
\keyword{phtotobacteria}
