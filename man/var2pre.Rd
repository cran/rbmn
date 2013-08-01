\name{var2pre}
\alias{var2pre}
\title{
returns the precision matrix from the variance
}
\description{
returns the precision matrix from the variance preserving possible
variable names
}
\usage{
var2pre(ma)
}
\examples{
 var2pre(rbmn0mn.04$gamma);
}
\arguments{
\item{ma}{ The variance matrix.}
}
\details{
Non full rank matrices are accepted, a generalized inverse is
returned and a warning is issued.
}
\value{
The precision matrix
}
