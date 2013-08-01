\name{mn4joint1condi}
\alias{mn4joint1condi}
\title{
computes a joint distribution from a marginal and a conditional one
for multinormal distributions
}
\description{
returns the expectation and variance of the multinormal normal
distribution defined through a marginal subcomponent and a
conditional distribution.
}
\usage{
mn4joint1condi(lmar, lcon)
}
\examples{
 lcon <- list(a=c(D=2, E=4), 
 b=matrix(1:6, 2, dimnames=list(LETTERS[4:5], 
 LETTERS[1:3])), 
 S=matrix(c(1, 1, 1, 2), 2));

 print8mn(mn4joint1condi(rbmn0mn.01, lcon));
}
\arguments{
\item{lmar}{ list defining the distribution of the marginal part with
\code{$mu}, its expectation, and \code{$gamma}, its variance matrix
(in fact a /mn/ object).}
\item{lcon}{ list defining the distribution of the conditional part
(see the \emph{Details} section).}
}
\details{
The conditional distribution is defined with a list having \code{$a}
for the constant part of the expectation; \code{$b} for the
regression coefficient part of the expectation; and \code{$S} for the
residual variance matrix.
}
\value{
A list: 

\item{\$mu}{The expectation vector.}
\item{\$gamma}{The joint variance matrix.}
 that is a /mn/ object.
}
