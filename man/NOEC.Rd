\name{NOEC}
\alias{NOEC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{NOEC and LOEC Calculation}
\description{calculating the NOEC and LOEC using Dunnett's test}
\usage{NOEC(x, rspn, blankC = FALSE, sigLev = 0.05, alternertive = 'B', sav = FALSE)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ a numeric vector of experimental concentrations}
  \item{rspn}{ a numeric matrix of experimental responses with at least three replicates.}
  \item{blankC}{TRUE if rspn contains responses of blank control. The default is FALSE.}
  \item{sigLev}{the significance level for Dunnett's test. The default is 0.05.}
  \item{alternertive}{the alternative hypothesis: "U"=upper one-sided test; "B"=two-sided 
	test(default).} 
  \item{sav}{TRUE: save output to a default file; FALSE: output will not be saved; 
    a custom file directory: save output to the custom file directory.} 
}
\details{ Dunnett's test (Dunnett, 1964) is performed to compare the treatment groups with the 
blank controls. The critical constants (store in DTcv) were calculated using step-down 
Dunnett test procedure. Three significance level (0.01, 0.05, and 0.1) are supported. 
	## Q: One dataset has four blank controls (C1, C2, C3, C4) and one treatment has three 
	replicates (T1, T2, T3), 
	##   another treatment has five replicates (R1, R2, R3, R4, R5), how to arrange the 
	response matrix (rspn)?
	## A: Label the missing values as NA, the response matrix (rspn) can be arranged as follows:\cr
		
		C1	C2	C3	C4	NA\cr
		T1	T2	T3 NA	NA\cr
	   	R1	R2	R3	R4	R5\cr
	 
	The adjustation of critical value for the unequal variances or unequal number of control 
	and replicates is skipped in this program. 
}
\value{
\item{mat }{information on Dunnett's test.DT: Dunnett's test values; DTcv: critical values 
for Dunnett's test at the significance level of sigLev.}
\item{noec }{non-observed effect concentration (NOEC).}
\item{loec }{least-observed effect concentration (LOEC).}
\item{sigLev }{the significance level used in the Dunnett's test.}
\item{DF }{the number of treatments and degree of fredom.}
\item{sav}{TRUE: save output to a default file; FALSE: output will not be saved; 
  a custom file directory: save output to the custom file directory.} 
}

\references{
Dunnett, C.W., 1964. New tables for multiple comparisons with a control. Biometrics 30, 482-491.
}

\note{x a vector of concentrations or levels in an ascending order.
	response matrix with at least 3 replicates.
	if the response matrix (rspn) contains blank controls (blankC = TRUE), 
	the blank controls should be allocated in the first tow of rspn matrix. 
	missing values should be labled as NA.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
## example 1
# calcualte the NOEC and LOEC of heavy metal Ni(2+) on the MCF-7 cells at the default significance 
# level of 0.05
x <- cytotox$Ni$x
rspn <- cytotox$Ni$y
NOEC(x, rspn)

## example 2
# calcualte the NOEC and LOEC of Neomycin sulfate on the phtotobacteria at the significance 
# level of 0.01
x <- antibiotox$NEO$x
rspn <- antibiotox$NEO$y
NOEC(x, rspn, sigLev = 0.01)
}
\keyword{NOEC}
\keyword{LOEC}
