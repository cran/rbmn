\name{print8nbn}
\alias{print8nbn}
\title{
print function for a /nbn/ object.
}
\description{
prints a /nbn/ object.
}
\usage{
print8nbn(nbn, what="pr", digits=3, ordering=NULL)
}
\examples{
 print8nbn(rbmn0nbn.01);
 print8nbn(rbmn0nbn.03, "pm", order=1:2)
}
\arguments{
\item{nbn}{\code{nbn} object to be printed.}
\item{what}{a \code{character(1)}; when comprising "p" the name of
each node with its parents are given, when comprising "r the formula
regression of each node is given with the node, when comprising "m"
the model is given.}
\item{digits}{ when not null, the number of digits for rounding.}
\item{ordering}{ Nodes are given following the indices of "ordering"
if \code{numeric} or the names if it is \code{character}. \code{NULL}
means the identity permutation. Repetitions or missing nodes are
accepted.}
}
\value{
Nothing but but \code{nbn} is printed.
}
