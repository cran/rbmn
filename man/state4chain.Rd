\name{state4chain}
\alias{state4chain}
\title{
returns the states of each node of a chain
}
\description{
From a \code{chain} object returns a named character precising the
role of each node: "r" for root, "c" for collider, "t" for
transmitter and "l" for leaf.
}
\usage{
state4chain(chain)
}
\examples{
 state4chain(rbmn0chain.01);
 state4chain(rbmn0chain.03);
}
\arguments{
\item{chain}{ the \code{chain} object to be considered.}
}
\value{
a character of the states named with node names.
}
