\name{marginal4chain}
\alias{marginal4chain}
\title{
returns marginal expectations and standard deviations of a chain
}
\description{
From a \code{chain} object returns a list with two components:
\code{$mu} and \code{$sigma} vectors of marginal expectations and
standard deviations.\cr
}
\usage{
marginal4chain(chain)
}
\examples{
 marginal4chain(rbmn0chain.02);
}
\arguments{
\item{chain}{ the \code{chain} object to be considered.}
}
\value{
a list with the two components \code{$mu} and \code{$sigma}.
}
