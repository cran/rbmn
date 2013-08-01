\name{generate8chain}
\alias{generate8chain}
\title{
generation of a /chain/ /nbn/
}
\description{
[randomly] generates a /chain/ /nbn/.
}
\usage{
generate8chain(rnn=c(3, 7), proo=0.5, rcor=c(-1, 1), rmu=c(0, 0), rsig=c(0, 1),
	  nona=r.form3names(max(rnn)))
}
\examples{
 set.seed(1234);
 print8chain(generate8chain());
 print8chain(generate8chain());
 print8chain(generate8chain(rnn=10, rcor=0.5));
 print8chain(generate8chain(rnn=10, rcor=0.5));
}
\arguments{
\item{rnn}{Range of the number of nodes.}
\item{proo}{Probabilit[y\eqn{|}ies] that the successive and
acceptable nodes be colliders. Can be a vector.}
\item{rcor}{Range of the correlations between neighbour nodes.}
\item{rmu}{Range of the expectations.}
\item{rsig}{Range of the standard deviations.}
\item{nona}{Proposed names for the maximum number of nodes, only the
necessary first ones will be used.}
}
\details{
Proposed ranges can be a unique value, implying no randomness in the
value.\cr Roots are placed according to \code{proo} probabilities,
then collider are placed in between with uniform probability on the
possibles nodes.
}
\value{
A /chain/ coding list is returned.
}
