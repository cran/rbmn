\name{adja2nbn}
\alias{adja2nbn}
\title{
standardized /nbn/ from an adjacency matrix
}
\description{
returns a \code{nbn} object with O/1 regression coefficients having
\code{adja} as adjacency matrix.
}
\usage{
adja2nbn(adja)
}
\examples{
 print8nbn(adja2nbn(adja4nbn(rbmn0nbn.03)));
}
\arguments{
\item{adja}{ The initial adjacency matrix.}
}
\value{
The corresponding standardized \code{nbn} object.
}
