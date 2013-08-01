\name{cor4var}
\alias{cor4var}
\title{
returns the correlation matrix from the variance
}
\description{
returns the correlation matrix from the variance preserving possible
variable names
}
\usage{
cor4var(ma)
}
\examples{
 cor4var(rbmn0mn.04$gamma);
}
\arguments{
\item{ma}{ The variance matrix.}
}
\details{
Zero variances are detected and accepted (all associated correlation
coefficients are forced to be zero.\eqn{>>}
}
\value{
The correlation matrix
}
