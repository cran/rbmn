\name{nbn2chain}
\alias{nbn2chain}
\title{
transforms a /nbn/ into a /chain/
}
\description{
returns the chain obtained from \code{nbn} which is supposed to a
chain. If it is not a chain, an error is issued.
}
\usage{
nbn2chain(nbn)
}
\examples{
 print8chain(nbn2chain(rbmn0nbn.02));
}
\arguments{
\item{nbn}{ The /nbn/ object to consider.}
}
\details{
It is advised to use \code{is8nbn8chain} before calling this
function.
}
\value{
The resulting chain
}
