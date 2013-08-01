\name{crossed4nbn1nbn}
\alias{crossed4nbn1nbn}
\title{
creates a crossed-nbn from two /nbn/s
}
\description{
A crossed /nbn/ is a /nbn/ obtained when replacing each node of the
first /nbn/ by the second /nbn/ and vice-versa.\cr Let \code{nn1/nn2}
and \code{na1/na2} be the node and arc numbers of the two
\code{nbn}s, the node number of the crossed \code{nbn} is
\code{nn1*nn2} and its arc number is \code{nn1*na2+nn2*na1}.\cr The
regression coefficients attributed to the crossed \code{nbn} are the
products of the weights (\code{we1/we2}) and the regression
coefficients of the initial \code{nbn}.
}
\usage{
crossed4nbn1nbn(nbn1, nbn2, we1=rep(1, length(nbn1)), we2=rep(1, length(nbn2)),
	  nona=as.vector(outer(names(nbn1), names(nbn2), paste,
	  sep="_")))
}
\examples{
 print8nbn(crossed4nbn1nbn(rbmn0nbn.01, rbmn0nbn.04));
}
\arguments{
\item{nbn1}{ The first generating /nbn/.}
\item{nbn2}{ The second generating /nbn/.}
\item{we1}{ The weight to apply to the nodes of the first generating
/nbn/.}
\item{we2}{ The weight to apply to the nodes of the second generating
/nbn/.}
\item{nona}{ The node names to give to the crossed /nbn/, the nodes
of the \code{nbn1} varying first.}
}
\details{
The \code{mu} coefficient is the sum of the two corresponding
\code{mu}s of the generating \code{nbn}.\cr The \code{sigma}
coefficient is the product of the two corresponding \code{sigma}s of
the generating \code{nbn}.\cr The regression coefficient are directed
inherited from the \code{nbn} which is duplicated with this arc.
}
\value{
The resulting crossed \code{nbn} object.
}
