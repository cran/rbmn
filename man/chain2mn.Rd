\name{chain2mn}
\alias{chain2mn}
\title{
computes the distribution of a chain
}
\description{
returns the /mn/ object associated to a /chain/ object. Much better
to use this function that the general function \code{nbn2mn} since
exact formulae are applied.
}
\usage{
chain2mn(chain, order=TRUE)
}
\examples{
 print8mn(chain2mn(rbmn0chain.01));
}
\arguments{
\item{chain}{ The chain object to consider.}
\item{order}{ Must a topological order be imposed?}
}
\value{
The resulting /mn/ object. Following the convention of \code{mn}
objects, a topological order is given to it. This is necessary to
retrieve the associate /nbn/.
}
