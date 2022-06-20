\name{DTcv}
\alias{DTcv}
\docType{data}
\title{critical value for Dunnett's test }
\description{ DTcv provides the critical constants calculated based step-down 
Dunnett test procedure. Three significance level (0.01, 0.05, and 0.1) each with two 
 alternative hypothesis ("U"=upper one-sided test; "B"=two-sided 
	test) are supported. 
}
\usage{data(DTcv)}
\format{ at most 30 treatments (1 : 30), 
		and 35 degree of freedom(c(5 : 30, 40, 50, 60, 80, 100, 120, 200, 1000, 3000)).
  \describe{
    \item{\code{DTcv}}{a matrix of critical value for Dunnett's test}
  }
}
\details{
> head(DTcv)\cr
\verb{   }df p twoside.01 twoside.05 twoside.10 oneside.01 oneside.05 oneside.10\cr
[1,] 5  1      4.032      2.571      2.015      3.365      2.015      1.476\cr
[2,] 5  2      4.627      3.030      2.433      3.900      2.440      1.873\cr
[3,] 5  3      4.948      3.294      2.669      4.225      2.681      2.095\cr
[4,] 5  4      5.218      3.474      2.831      4.434      2.848      2.245\cr
[5,] 5  5      5.416      3.616      2.956      4.585      2.974      2.360\cr
[6,] 5  6      5.538      3.727      3.055      4.723      3.080      2.451\cr
df: degree of freedom; p: the number of treatment.
}
\source{ The critical constants (store in DTcv) were calculated using step-down 
Dunnett test procedure( the cvSDDT funtion in R package DunnettTests).}

\keyword{critical value}
\keyword{degree of freedom}
\keyword{significance level}
