\name{gema2mn}
\alias{gema2mn}
\title{
computes a /mn/ from a /gema/
}
\description{
from a /gema/ object defining a normal Bayesian network, computes the
expectation and variance matrix.
}
\usage{
gema2mn(gema)
}
\examples{
 print8mn(gema2mn(rbmn0gema.04));
}
\arguments{
\item{gema}{Initial \code{gema} object.}
}
\value{
a list with the following components: \code{mu} and \code{gamma}.
}
