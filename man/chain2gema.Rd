\name{chain2gema}
\alias{chain2gema}
\title{
transforms a /chain/ to a /gema/
}
\description{
From a \code{chain} object returns the \code{gema} using a direct
formulae.\cr Much precised than to use the /nbn/ way.
}
\usage{
chain2gema(chain)
}
\examples{
 identical(chain2gema(rbmn0chain.02)$mu, rbmn0gema.02$mu);
 print(chain2gema(rbmn0chain.02)$li-rbmn0gema.02$li);
}
\arguments{
\item{chain}{ the \code{chain} object to be transformed.}
}
\value{
The corresponding \code{gema} object.
}
