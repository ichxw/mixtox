\name{figPlot}
\alias{figPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Ploting concentration response curve}
\description{Plotting concentration response curves with/without confidence intervals.
}
\usage{figPlot(object, logT = TRUE, xlabel = "concentration (mol/L)", ylabel = "Response", 
               ylimit, lgd = NULL)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{object of class curveFit.}
  \item{logT}{logarithm transformation on concentration: TRUE or FALSE(defalt if TRUE).}
  \item{xlabel, ylabel}{plot labels.}
  \item{ylimit}{range of the y axis.}
  \item{lgd}{legend of the plot.}
}
\details{Plot the concentration response curves, experimental data, CI, and PI in one figure.
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\examples{
## example 1
# 
x <- antibiotox$PAR$x
expr <- antibiotox$PAR$y
obj <- curveFit(x, expr, eq = 'Logit', rtype = 'quantal', param = c(26, 4), effv = c(0.05, 0.5))
figPlot(obj)

## example 2
# 
x <- hormesis$HmimCl$x
rspn <- hormesis$HmimCl$y
obj <- curveFit(x, rspn, eq = 'Biphasic', param = c(-0.59, 0.001, 160,0.05, 19),
         effv = c(0.05, 0.5), rtype = 'hormesis')
figPlot(obj, logT = TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{concentration response curve}
