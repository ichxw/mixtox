\name{readModel}
\alias{readModel}

\title{Read curve fitting information  }
\description{ Read curve fitting information in specific format from a txt or csv file
}
\usage{readModel(File)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{File}{a txt or csv file with curve fitting information. a txt file should be tab or space separated; 
       a csv file should be comma separated.}


}

\format{ a file contains curve fitting information\cr 
The file extension needs to be .txt if it's a space or tab separated file.\cr
The file extension needs to be .csv if it's a comma-separated values file.\cr
}
\details{The file needs three paramters (alpha	beta	gamma), just set gamma to 0 for 
         models with only 2 parameters(e.g., Logit/Weibull)\cr
> The format of the file readModel can read:\cr
ID	Model	alpha	beta	gamma\cr
Cu	Logit	9.91	22.15	0\cr
PAR	Weibull	18.31	3.32	0\cr

}

\keyword{curve fitting}
\keyword{txt}
\keyword{csv}
