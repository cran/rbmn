\name{check8gema}
\alias{check8gema}
\title{
checks a /gema/ object
}
\description{
checks the consistency of \code{gema} as a /gema/ object issues a
fatal error with some clues if inconsistent.
}
\usage{
check8gema(gema)
}
\examples{
 check8gema(rbmn0gema.01);
 res <- check8gema(rbmn0adja.01);
 if (is.na(as.logical(res))) { print(res);}
}
\arguments{
\item{gema}{ The \code{gema} object to check.}
}
\details{
Looking a the code of this function provides a way to know which are
the requirements of a /chain/ object.
}
\value{
\code{TRUE} or a \code{character} containing some clue about the
discovered inconsistency.
}
