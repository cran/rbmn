\name{arcs4nbn1nbn}
\alias{arcs4nbn1nbn}
\title{
returns the list of 'parallel' arcs of a crossed-nbn
}
\description{
Returns a list of matrices with two columns (as needed by
\code{estimate8constrainednbn}) indicating corresponding arcs for
each arcs/nodes of \code{nbn1} (or \code{nbn2}) of the crossed /nbn/
obtained when crossing /nbn1/ and /nbn2/ with node names given by
\code{nona}.
}
\usage{
arcs4nbn1nbn(nbn1, nbn2, type="a1", nona=as.vector(outer(names(nbn1),
	  names(nbn2), paste, sep="_")))
}
\examples{
 print(arcs4nbn1nbn(rbmn0nbn.01, rbmn0nbn.04));
}
\arguments{
\item{nbn1}{ The first generating /nbn/.}
\item{nbn2}{ The second generating /nbn/.}
\item{type}{ Must be \code{"a1"} to indicate that the parallelism
must be done for each arc of \code{nbn1}. \code{"a2"} for each arc of
\code{nbn2}. Or \code{"n1"} for each node of \code{nbn1}. Or
\code{"n2"} for each node of \code{nbn2}. }
\item{nona}{ The node names to give to the crossed /nbn/, the nodes
of the \code{nbn1} varying first.}
}
\value{
The resulting named (after node names) list of matrices.
}
