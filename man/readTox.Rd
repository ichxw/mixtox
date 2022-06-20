\name{readTox}
\alias{readTox}

\title{Read dose-response data  }
\description{ Read dose-response data of specific format from a txt or csv file
}
\usage{readTox(File, light=TRUE)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{File}{a txt or csv file with dose-response data. a txt file should be tab or space separated; 
       a csv file should be comma separated.}
  \item{light}{TRUE: the file only contrains doses and responses (light version); 
    FALSE: the file also contains information about experiment design and data type 
    (complex version)}

}

\format{ Dose-response data contains the following information: the name of compound,
the number of doses, the repititions, and the data type.\cr 
The file extension needs to be .txt if it's a space or tab separated file.\cr
The file extension needs to be .csv if it's a comma-separated values file.\cr
}
\details{
> The light format of the file readTox can read:\cr
conc	tier1	tier2	tier3\cr
2.50E-07	0.06	0.03	0.07\cr
3.50E-07	0.04	0.07	0.04\cr
5.34E-07	0.07	0.15	0.06\cr
7.68E-07	0.16	0.14	0.10\cr
1.13E-06	0.21	0.26	0.23\cr
1.67E-06	0.34	0.30	0.31\cr
2.50E-06	0.46	0.44	0.48\cr
3.50E-06	0.52	0.59	0.61\cr
5.34E-06	0.73	0.69	0.71\cr
7.68E-06	0.79	0.67	0.78\cr
1.13E-05	0.80	0.72	0.81\cr
1.67E-05	0.82	0.80	0.82\cr

> The complex format of the file readTox can read:\cr
cu	12	3	type\cr	
ID	conc	tier1	tier2	tier3\cr
1	2.50E-07	0.06	0.03	0.07\cr
2	3.50E-07	0.04	0.07	0.04\cr
3	5.34E-07	0.07	0.15	0.06\cr
4	7.68E-07	0.16	0.14	0.10\cr
5	1.13E-06	0.21	0.26	0.23\cr
6	1.67E-06	0.34	0.30	0.31\cr
7	2.50E-06	0.46	0.44	0.48\cr
8	3.50E-06	0.52	0.59	0.61\cr
9	5.34E-06	0.73	0.69	0.71\cr
10	7.68E-06	0.79	0.67	0.78\cr
11	1.13E-05	0.80	0.72	0.81\cr
12	1.67E-05	0.82	0.80	0.82\cr
}

\keyword{dose-response data}
\keyword{txt}
\keyword{csv}
