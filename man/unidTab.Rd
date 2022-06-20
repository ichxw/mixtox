\name{unidTab}
\alias{unidTab}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Uniform Design Table}
\description{The good lattice point method with a power generator 
			was used to construct the uniform 
			experimental tables . The centered L2-discrepancy (CD2) and the symmetric discrepancy 
			algorithms (sd2) were employed to measure the uniformity and 
			find the one with lowest discrepancy.
			}
\usage{unidTab(lev, fac, algo = "cd2", sav = FALSE)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lev}{the number of runs (levels or pseudo-levels)}
  \item{fac}{the number of factors}
  \item{algo}{algorithms used to calculate the discrepancy. "cd2": the centered L2-discrepancy 
		algorithm. "sd2": the symmetric discrepancy algorithm.
		}
  \item{sav}{TRUE: save output to a default file; FALSE: output will not be saved; 
  a custom file directory: save output to the custom file directory.} 
}
\details{Uniform design (UD) was proposed by Fang et al (Acta Math Appl Sin 3:363-372 (1980)). 
		An appropriate uniform design table is constructed 
		according to the factor (the number mixture 
		components) and level (the number of experiments need to run). 
		Many methods can be used to 
		construct the uniform table. In the past decades many methods have been proposed for 
		constructing (nearly) uniform designs, such as the good lattice point (glp) method, the 
		glp method with a power generator (pglp method) (Fang 1980; Fang andWang 1994), 
		the cutting method (Maand Fang 2004), the optimization method (Winker and Fang 1998).\cr
		However, when s is large, the glp method has a large computational cost. And the pglp 
		method has the lowest computation complexity among various methods in quasi-Monte Carlo 
		methods and a good performance when \eqn{n} or \eqn{n + 1} is a prime number 
		and \eqn{s} is small (Fang 1980; 
		Fang and Wang 1994), while the pglp method may have a poor performance when s is large.
		Here, we choose the glp method with a power generator to construct the uniform table. The 
		centered L2-discrepancy (cd2) is set as default over the symmetric discrepancy algorithm
		for its accuracy. The cd2 algorithm is defined as follows:
		\deqn{C{D_2}(P) = {\left[ {{{\left( {\frac{{13}}{{12}}} \right)}^s} - \frac{{{2^{1 - s}}}}
		{n}\sum\limits_{k = 1}^n {\prod\limits_{i = 1}^s {{\theta _{ki}} + \frac{1}{{{n^2}}}
		\sum\limits_{k,l = 1}^n {\prod\limits_{i = 1}^s {{\phi _{k,li}}} } } } } 
		\right]^{\frac{1}{2}}}}
		With the definition of \eqn{\theta _{ki}} and \eqn{\phi _{k,li}} as follows:
		\deqn{{\theta _{ki}} = 2 + \left| {{x_{ki}} - \frac{1}{2}} \right| - 
		{\left| {{x_{ki}} - \frac{1}{2}} \right|^2}}
		\deqn{{\phi _{k,li}} = 1 + \frac{1}{2}\left( {\left| {{x_{ki}} - \frac{1}{2}}
		\right| + \left| {{x_{li}} - \frac{1}{2}} \right| - \left| {{x_{ki}} - 
		{x_{li}}} \right|} \right)}
		where \eqn{n}, \eqn{s} are the number of runs (levels or multiple of levels) and 
		the number of input variables (factors), respectively.		
		}
\value{
  \item{T }{the selected uniform table(s) . It may contain one or more uniform 
		tables with the same discrepancy
		}
  \item{D }{the discrepancy of the constructed uniform table}
  \item{sav }{TRUE: save output to a default file; FALSE: output will not be saved; 
  a custom file directory: save output to the custom file directory.} 
}
\references{
Wang, Y., Fang, K.-T., 1996. Uniform design of experiments with mixtures. Sci. China Ser. 
A-Mathematics Phys. Astron. 39, 264-275.\cr
Hickernell, F.J., 1996. A generalized discrepancy and quadrature error bound. Math. Comput. 
67, 299-322.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~
\examples{
## example 1
# construct uniform table with 11 runs and 7 factors using the default centered L2-discrepancy 
# algorithm
unidTab(11, 7)

## example 2
# construct uniform table with 37 runs and 13 factors using the symmetric discrepancy algorithm 
unidTab(lev = 37, fac = 13, algo = "sd2")

## example 3
# construct uniform table with 37 runs and 13 factors using default centered L2-discrepancy 
# algorithm
unidTab(lev = 37, fac = 13, algo = "cd2")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{uniform design}
\keyword{ good lattice point }
