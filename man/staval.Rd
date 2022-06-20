\name{staval}
\alias{staval}
\docType{data}
\title{Starting Values for 13 Sigmoidal and  4 Hormetic Models}
\description{providing staring values for function tuneFit. \cr
}
\usage{data(staval)}
\format{
  A list with starting values for 13 sigmoidal models and six hormetic models
  \describe{
    \item{\code{staval$Hill}}{sarting values for Hill model}
    \item{\code{staval$GL}}{sarting values for generalized logit model}
	\item{\code{staval}}{starting values for 13 sigmoidal and  five hormetic models }	
  }
}
\details{
Thirteen monotonic(sigmoidal) models ("Hill", "Hill_two", "Hill_three", "Hill_four",
"Weibull", "Weibull_three", "Weibull_four", "Logit", "Logit_three", "Logit_four",
"BCW(Box-Cox-Weibull)", "BCL(Box-Cox-Logit)", "GL(Generalized Logit)") and four
non-monotonic(J-shaped) models ("Brain_Consens", "BCV", "Biphasic", "Hill_five")
}
\examples{
# example 1
## Retrive the starting values for Hill.
staval$Hill

# example 2
## Retrive the starting values for Weibull.
staval$Weibull
}
\keyword{starting values}
\keyword{curve fitting}
