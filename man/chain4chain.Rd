\name{chain4chain}
\alias{chain4chain}
\title{
extracts a chain from a chain
}
\description{
returns the chain obtained from \code{chain} retaining only nodes
indicated by \code{nodes} and conditioned with nodes indicated in
\code{condi}.
}
\usage{
chain4chain(chain, nodes, condi=numeric(0), value=rep(0, length(condi)))
}
\examples{
 chain4chain(rbmn0chain.02, c("a", "d"), c("b"), 12);
}
\arguments{
\item{chain}{ The chain object to consider.}
\item{nodes}{ \code{numeric} (or \code{character}) vector giving the
numbers (or names) of the nodes to be retained in the extracted
chain.}
\item{condi}{ \code{numeric} (or \code{character}) vector giving the
numbers (or names) of the conditioning nodes for the extracted
chain.}
\item{value}{ Numerical values associated to \code{condi}.}
}
\details{
Integration is done for nodes not belonging to the extracted chain
nor being in the conditioning subset. Then the distribution of the
retained nodes is left identical to this in the initial chain.
}
\value{
The resulting chain
}
