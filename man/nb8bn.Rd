\name{nb8bn}
\alias{nb8bn}
\title{
number of Bayesian networks
}
\description{
returns the number of different Bayesian networks having \code{n}
labelled or not nodes. Non labelled nodes means that nodes are
exchangeable: \code{A -> B} is identical to \code{A <- B}.
}
\usage{
nb8bn(n, label=FALSE)
}
\examples{
 nb8bn(5)
 nb8bn(5, TRUE);
}
\arguments{
\item{n}{ number of nodes. Must be less or equal to 18.}
\item{label}{Indicates if the nodes must be considered as labelled or
not.}
}
\details{
When not labelled nodes, the results were proposed by Sloane in 'the
on line encyclopedy of integer sequences' (http://oeis.org/A003087).
For labelled nodes, just the application of the recursive formula of
Robinson.
}
\value{
Number of Bayesian networks
}
