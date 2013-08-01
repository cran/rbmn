\name{rm8nd4nbn}
\alias{rm8nd4nbn}
\title{
removes some nodes from a /nbn/
}
\description{
returns a /nbn/ object deduced from an original /nbn/ by integrating
on a given subset of nodes.
}
\usage{
rm8nd4nbn(nbn, nodes)
}
\examples{
 rm8nd4nbn(rbmn0nbn.04, "1.1"); 
}
\arguments{
\item{nbn}{ The \code{nbn} object to reduce.}
\item{nodes}{\code{character} or \code{numeric} vector giving the
subset of nodes to remove.}
}
\details{
The transformation is made through the associated joint distributions
for the probabilities and with the help of the function
\code{rm8nd4adja} for the relationships.
}
\value{
The resulting \code{nbn}.
}
