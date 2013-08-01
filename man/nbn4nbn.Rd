\name{nbn4nbn}
\alias{nbn4nbn}
\title{
From a /nbn/ computes the associated nbn1
}
\description{
returns a /nbn/ object with the same structure as \code{nbn} but all
\code{$mu} are put to zero, all \code{$sigma} to one as well as
\code{$regcof}.
}
\usage{
nbn4nbn(nbn)
}
\examples{
 print8nbn(nbn4nbn(rbmn0nbn.04));
}
\arguments{
\item{nbn}{ The \code{nbn} object to transform.}
}
\details{
These coefficient values allows the easy study of the /nbn/
structure.
}
\value{
The resulting \code{nbn}.
}
