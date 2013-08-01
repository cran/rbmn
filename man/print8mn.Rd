\name{print8mn}
\alias{print8mn}
\title{
standard print function for a /mn/ object.
}
\description{
prints a /mn/ object completely or a part of it.
}
\usage{
print8mn(mn, what="msC", ordering=NULL, digits=3, printed=TRUE)
}
\examples{
 print8mn(rbmn0mn.01);
}
\arguments{
\item{mn}{\code{mn} object to be printed.}
\item{what}{a \code{character(1)}; when comprising "m" the
expectations are printed, "s" the standard deviations are printed,
"C" the correlation matrix is printed, "S" the variance matrix is
printed, "P" the precision matrix is printed, "p" the normalized
precision matrix is printed.}
\item{ordering}{ Nodes are given following the indices of "ordering"
if \code{numeric} or the names if it is \code{character}. \code{NULL}
means the identity permutation. Repetitions or missing nodes are
accepted.}
\item{digits}{ when not null, the number of digits for rounding the
parameter values.}
\item{printed}{ \code{TRUE} to issue a printing, if not the prepared
matrix is returned.}
}
\value{
The \code{mn} is printed or a matrix having \code{nn x ?} is returned
binding which elements precised in the argument \code{what}.
}
