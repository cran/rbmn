\name{estimate8nbn}
\alias{estimate8nbn}
\title{
estimating the /nbn/ parameters
}
\description{
From a /nbn/ to describe the DAG, and a data.frame containing the
necessary observations, returns the /nbn/ with all its parameters
newly estimated.
}
\usage{
estimate8nbn(nbn, data)
}
\examples{
 data(boco);
 print8nbn(rbmn0nbn.05);
 print8nbn(estimate8nbn(rbmn0nbn.05, boco));
}
\arguments{
\item{nbn}{ The initial /nbn/.}
\item{data}{The data frame comprising all /nbn/ nodes.}
}
\details{
No constraints are put on the parameters.
}
\value{
The resulting /nbn/ with the estimated parameters.
}
