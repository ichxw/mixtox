\name{qq4res}
\alias{qq4res}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Residual Normal QQ Plot }
\description{Producing a side-by-side QQ plot of the residuals against standard normal quantiles.
}
\usage{qq4res(object, xlabel = 'Theoretical Quantiles', ylabel = 'Residuals', 
       lgd = NULL)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{object of class curveFit.}
  \item{xlabel, ylabel}{plot labels.}
  \item{lgd}{legend of the plot}
}
\details{The empirical quantiles are plotted against the quantiles of a standard normal 
		distribution. If the residuals are from a normal distribution with mean 0, the points 
		tend to fall along the reference line that has an intercept of 0 and a slope equal to 
		the estimated standard deviation.
}


%% ~Make other sections like Warning with \section{Warning }{....} ~
\examples{
## example 1
# 
x <- antibiotox$PAR$x
expr <- antibiotox$PAR$y
obj <- curveFit(x, expr, eq = 'Logit', rtype = 'quantal', param = c(26, 4), effv = c(0.05, 0.5))
qq4res(obj)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{normal QQ plot}
