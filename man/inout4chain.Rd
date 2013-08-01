\name{inout4chain}
\alias{inout4chain}
\title{
reduces a chain to its inputs and outputs
}
\description{
From a \code{chain} returns the reduced \code{chain} comprising only
inputs (that is root nodes) and outputs (that is colliders and ends
which are not roots)
}
\usage{
inout4chain(chain)
}
\examples{
 print8chain(inout4chain(rbmn0chain.02));
}
\arguments{
\item{chain}{ The chain object to consider.}
}
\value{
The resulting chain
}
