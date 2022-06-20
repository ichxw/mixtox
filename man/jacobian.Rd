\name{jacobian}
\alias{jacobian}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Jacobian Matrix Calculation}
\description{calculating the Jacobian matrix for confidence intervals}
\usage{jacobian(eq, x, paraHat)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eq}{ equation}
  \item{x}{ numeric vector of experimental concentrations}
  \item{paraHat}{ fitted parameters}
 
}
\details{ 
Jacobian matrix is the matrix of all first-order partial derivatives of a vector-valued function.
}
\value{
\item{jac }{Jacobian matrix.}

}
\references{ https://en.wikipedia.org/wiki/Jacobian_matrix_and_determinant
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
## example 1
x <- cytotox$Ni$x
rspn <- cytotox$Ni$y
obj <- curveFit(x, rspn, eq = 'Logit', param = c(12, 3), effv = c(0.05, 0.5), rtype = 'quantal')
jacobian('Logit', x, obj$p)

}
\keyword{Jacobian matrix}

