\name{hormesis}
\alias{hormesis}
\docType{data}
\title{Non-monotonic Concentration-response Data}
\description{
Two ionic liquids: 1-Octyl-3-methylimidazolium chloride (Omim) and 
1-Hexyl-3-Methylimidazolium tetrafluoroborate (Hmim). 
Two organic solvents: Acetonitrile(ACN) and acetone.\cr
The concentration-response data include the effect of these four compounds 
on firefly luciferase. 
Those concentration-responses were fitted using non-monotonic models. 
}
\usage{data(hormesis)}
\format{
  A list with non-monotonic concentration-responses of four chemicals and associated 
  fitting information.
  \describe{
    \item{\code{OmimCl$x}}{a numeric vector of concentrations}
    \item{\code{OmimCl$y}}{a numeric matrix of responses}
	\item{\code{OmimCl$name}}{name of test substance}
	\item{\code{OmimCl$concNum}}{the number of concentrations}
	\item{\code{OmimCl$tierNum}}{the number of repetitions}
	\item{\code{OmimCl$type}}{type of test substance}
	
	\item{\code{sgl$model}}{model used to fit the concentration-response data}
	\item{\code{sgl$param}}{fitted coefficients of those curves corresponding 
		to sgl$model}
	\item{\code{sgl$minx}}{numeric vector with  multiple concentrations that 
               induce maximun stimulation}
	\item{\code{sgl$miny}}{the largest stimulation}	
  }
}
\details{
The non-monotonic concentration-responses need to be scaled into [0, 1] using the 
following equation:
\deqn{E = \frac{{{I_0} - {I_i}}}{{{I_0}}}}
where \eqn{I_0} is the average of controls for inhibition test or the average 
of the maximum effect for the stimulation test and \eqn{I_i} the average 
effect of the \eqn{i^{th}} treatment.
}
\source{ Eexperiments were conducted in our lab. Detailed descripttion can be found in references }
\references{
Zhu X-W, et.al. 2013. Modeling non-monotonic dose-response relationships: Model evaluation
and hormetic quantities exploration. Ecotoxicol. Environ. Saf. 89:130-136.
}
\examples{
#example 1
## Retrive the toxicity data of acetonitrile on firefly luciferase.
hormesis$ACN

#example 2
## Retrive the minx of OmimCl, HmimCl, ACN, and Acetone
hormesis$sgl$minx
}
\keyword{hormesis}
\keyword{non-monotonic concentration-response data}