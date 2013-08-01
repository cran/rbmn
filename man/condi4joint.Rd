\name{condi4joint}
\alias{condi4joint}
\title{
computes some conditional distribution of a multinormal vector
}
\description{
returns the expectation and variance of a sub-vector conditionned
with another (non overlapping) sub-vector from an initial random
vector described by \code{mn}.
}
\usage{
condi4joint(mn, par, pour, x2=rep(0, length(pour)))
}
\examples{
 print8mn(condi4joint(rbmn0mn.04, c("1.1", "2.2", "1.2", "2.1"), NULL));
 print8mn(condi4joint(rbmn0mn.04, c("1.1", "2.2", "1.2", "2.1"), "C", 0));
 print(condi4joint(rbmn0mn.04, c("1.1", "2.2", "1.2", "2.1"), "C", NULL));
}
\arguments{
\item{mn}{ list defining the distribution of the initial vector with
\code{$mu}, its expectation, and \code{$gamma}, its variance matrix.}
\item{par}{ names (or indices) of the sub-vector to give the
distribution.}
\item{pour}{ names (or indices) of the conditionning sub-vector (can
be \code{NULL} when for non conditionning.}
\item{x2}{values to consider for the conditioning sub-vector. When
\code{NULL} the general form is supplied, not a /mn/ object. }
}
\details{
when no names are given to \code{mn$mu}, \code{par} and \code{pour}
are supposed containing indices and default sequential names are
provided.
}
\value{
A list:\cr when \code{x2} provides the values taken by the
conditioning part, it is a /mn/ object with its two components:
\code{$mu} for the expectation vector and \code{$gamma} for the
variance matrix.\cr when \code{x2} is \code{NULL} the list has got
three components: \code{$a} for the fixed part of the expectation
vector, \code{$b} for the regression coefficients to be associated to
the non precised \code{x2} values, varying part of the expectation
and \code{$gamma} for the variance matrix.\cr
}
