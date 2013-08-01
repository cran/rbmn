\name{nbn2nbn}
\alias{nbn2nbn}
\title{
computes the /nbn/ changing its topological order
}
\description{
returns the proposed /nbn/ with a new topological order without
modifying the joint distribution of all variables.\cr This allows to
directly find regression formulae within the Gaussian Bayesian
networks.
}
\usage{
nbn2nbn(nbn, norder)
}
\examples{
 print8mn(nbn2mn(rbmn0nbn.01, algo=1));
 print8mn(nbn2mn(rbmn0nbn.01, algo=2));
 print8mn(nbn2mn(rbmn0nbn.01, algo=3));
 print8mn(nbn2mn(nbn2nbn(rbmn0nbn.02, c(1, 2, 4, 5, 3))));
 print8mn(nbn2mn(nbn2nbn(rbmn0nbn.02, c(4, 1, 2, 3, 5))));
}
\arguments{
\item{nbn}{ The /nbn/ to transform.}
\item{norder}{ The topological order to follow. It can be indicated
by names or numbers. When not all nodes are included, the resulting
/nbn/ is restricted to these nodes after marginalization.}
}
\details{
BE aware that for the moment, no check is made about the topological
order and if it is not, the result is FALSE!
}
\value{
The resulting /nbn/.
}
