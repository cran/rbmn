\name{generate8nbn}
\alias{generate8nbn}
\title{
returns a randomly built /nbn/ object.
}
\description{
To obtain systematic results, you have to call \code{set.seed} before
hands.
}
\usage{
generate8nbn(rnn=c(3, 7), ppar=0.5, rreg=c(-1, 1), rmu=c(0, 0), rsig=c(0, 1),
	  nona=r.form3names(max(rnn)))
}
\examples{
 set.seed(1234)
 print8nbn(generate8nbn());
 print8nbn(generate8nbn());
}
\arguments{
\item{rnn}{Range of the number of nodes.}
\item{ppar}{Probabilities (not a range) of the parent occurrence for
each ancestor of every node. Can be a vector, cycled as necessary.}
\item{rreg}{Range of regression coefficients.}
\item{rmu}{Range of the conditional expectations.}
\item{rsig}{Range of the conditional standard deviations.}
\item{nona}{Proposed names for the maximum number of nodes, only the
necessary first ones will be used.}
}
\details{
Node numbers are uniformly drawn. Parent numbers are independently
drawn from all ancestors with the probability associated to the
considered node. Regression coefficient are uniformly drawn.
Conditional expectations and standard deviations are uniformly
drawn.\cr All range arguments can be given one value instead of two,
to precise the unique value to use.
}
\value{
a /nbn/ object, with nodes in topological order.
}
