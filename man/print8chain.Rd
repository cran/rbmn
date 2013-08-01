\name{print8chain}
\alias{print8chain}
\title{
prints a /chain/ object
}
\description{
prints a /chain/ object.
}
\usage{
print8chain(chain, digits=3)
}
\examples{
 print8chain(rbmn0chain.01);
 print8chain(rbmn0chain.02);
 print8chain(rbmn0chain.03);
}
\arguments{
\item{chain}{ The \code{chain} object to print.}
\item{digits}{ when not null, the number of digits for rounding the
numerical values.}
}
\details{
See \code{nbn2chain} code for some details about the definition of a
/chain/.
}
\value{
nothing but something is printed
}
