\name{nbn2rr}
\alias{nbn2rr}
\title{
computes standard matrices from a /nbn/
}
\description{
from a /nbn/ object defining a normal Bayesian network, returns a
list comprising (i) \code{mm} the vector of the mean of the different
nodes when the parents are nought, (ii) \code{ss} the vector of the
conditional standard deviations and (iii) \code{rr} the matrix of the
regression coefficients of the direct parents (\code{rr[i,j]}
contains the regression coefficient of the node \code{j} for its
parents \code{i} or zero when \code{i} is not a parent of \code{j}.
}
\usage{
nbn2rr(nbn)
}
\examples{
 nbn2rr(rbmn0nbn.01);
}
\arguments{
\item{nbn}{\code{nbn} object.}
}
\value{
the resulting list with the three components: \code{mm}, \code{ss}
and \code{rr}.
}
