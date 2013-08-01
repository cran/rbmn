\name{print8gema}
\alias{print8gema}
\title{
standard print function for a /gema/ object.
}
\description{
prints a /gema/ object completely or a part of it according to
\code{what} specification.
}
\usage{
print8gema(gema, what="ml", ordering=NULL, digits=3, printed=TRUE)
}
\examples{
 print8gema(rbmn0gema.01);
 print8gema(rbmn0gema.02, "m");
 print8gema(rbmn0gema.03, "l", digit=1);
 print8gema(rbmn0gema.04, printed=FALSE);
}
\arguments{
\item{gema}{\code{gema} object to be printed.}
\item{what}{a \code{character(1)}; when comprising "m" the
expectations are printed, "l" the linear combinations are printed.}
\item{ordering}{ Nodes are given following the indices of "ordering"
if \code{numeric} or the names if it is \code{character}. \code{NULL}
means the identity permutation. Repetitions or missing nodes are
accepted.}
\item{digits}{ when not null, the number of digits for rounding.}
\item{printed}{ \code{TRUE} to issue a printing, if not the prepared
matrix is returned.}
}
\value{
The \code{gema} is printed or a matrix having \code{nn x ?} is
returned binding which elements are precised in the argument
\code{what}.
}
