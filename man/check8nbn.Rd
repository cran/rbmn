\name{check8nbn}
\alias{check8nbn}
\title{
checks a /nbn/ object
}
\description{
checks the consistency of \code{nbn} as a /nbn/ object issues a fatal
error with some clues if inconsistent.
}
\usage{
check8nbn(nbn)
}
\examples{
 check8nbn(rbmn0nbn.01);
 res <- check8nbn(rbmn0adja.01);
 if (is.na(as.logical(res))) { print(res);}
}
\arguments{
\item{nbn}{ The \code{nbn} object to check.}
}
\details{
Looking a the code of this function provides a way to know which are
the requirements of a /chain/ object.
}
\value{
\code{TRUE} or a \code{character} containing some clue about the
discovered inconsistency.
}
