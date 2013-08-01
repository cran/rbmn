\name{mn2gema}
\alias{mn2gema}
\title{
computes a /gema/ from a /mn/
}
\description{
proposes generating matrices of a Bayesian network from a /mn/ object
defining a multinormal distribution by expectation and variance,
under the assumption that the nodes are in topological order.
}
\usage{
mn2gema(mn)
}
\examples{
 print8gema(mn2gema(rbmn0mn.04));
}
\arguments{
\item{mn}{Initial \code{mn} object.}
}
\value{
a list with the /gema/ components \code{$mu} and \code{$li}.
}
