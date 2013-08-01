\name{chain2correlation}
\alias{chain2correlation}
\title{
computes the correlation matrix of a chain
}
\description{
returns the correlation matrix of a /chain/ object.
}
\usage{
chain2correlation(chain)
}
\examples{
 chain2correlation(rbmn0chain.03);
}
\arguments{
\item{chain}{ The chain object to consider.}
}
\value{
The correlation matrix. It is not sorted to respect a topological
order contrary to \code{chain2mn} function.
}
