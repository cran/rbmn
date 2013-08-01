\name{gema2nbn}
\alias{gema2nbn}
\title{
computes a /nbn/ from a /gema/
}
\description{
from a /gema/ object defining a normal Bayesian network, computes
more standard /nbn/ where each node is defined from its parents.
}
\usage{
gema2nbn(gema)
}
\examples{
 print8nbn(gema2nbn(rbmn0gema.02));
}
\arguments{
\item{gema}{Initial \code{gema} object.}
}
\details{
using general formulae rather a sequential algorithm as done in the
original \code{gema2nbn} implementation.
}
\value{
the corresponding /nbn/.
}
