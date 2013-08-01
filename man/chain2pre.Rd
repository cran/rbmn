\name{chain2pre}
\alias{chain2pre}
\title{
computes the precision of a chain
}
\description{
returns the precision matrix of a chain, that is the inverse of its
variance (correlation) matrix. Much better to use this function that
\code{solve(chain2mn(chain)$gamma)} since exact formulae are applied.
}
\usage{
chain2pre(chain, corre=FALSE)
}
\examples{
 chain2pre(rbmn0chain.02);
}
\arguments{
\item{chain}{ The chain object to consider.}
\item{corre}{To get the inverse of the correlation matrix instead
of.}
}
\value{
A dimnamed matrix
}
