\name{rbmn-package}
\alias{rbmn-package}
\alias{rbmn}
\docType{package}
\title{
  Linear Gaussian Bayesian network manipulations
}

\description{
  General functions to generate, transform, display general and
  particular linear Gaussian Bayesian networks [/nbn/] are provided.\cr
  Specific /nbn/ are chain and crossed /nbn/s. Focus is given in getting
  joint and conditional probability distributions of the set of
  nodes.\cr
  \emph{rbmn} stands for \emph{R'eseau Bay'esien MultiNormal}.
}

\details{

  Some basic concepts:
  
  \itemize{
    \item \emph{chain /nbn/s} are /nbn/s where all nodes are connected with two
    other nodes, except the two ending nodes of the chain having only one
    connection. (This is not the usual terminology in graphical models
    but I didn't find a more appropriate word: suggestions are welcome.)
    \item \emph{crossed /nbn/s} are /nbn/s having the node set defined as a
    Cartesian product of two series of items, and a DAG based on this
    structure. See the \code{crossed4nbn1nbn} function and/or Tian (2013) for
    details.
    \item An \emph{adjacency matrix} is a matrix equivalent to the DAG
    associated to a /nbn/. Its rows as well as its columns are associated
    to the set of nodes. The \code{(i,j)} cell is one when there is an arc going
    from node \code{i} to \code{j} and zero otherwise.
  }

  Three equivalent ways can be used to represent the joint probability
  distribution of a set of nodes respectively associated to the
  structures /mn/, /nbn/ and /gema/:

  \itemize{

    \item /mn/ (for multivariate normal) is just the list of the
    expectation (\code{$mu}) and the variance matrix (\code{$gamma}).
    \item /nbn/ (for normal Bayesian network) is a simple list, a
    component a node described with a list. The names are node names and
    each list associated to a node provides the conditional expectation
    and variance, the parent (if any) and the associated regression
    coefficients.
    \item /gema/ (for generating matrices) is a list of a vector
    (\code{mu}) and a matrix (\code{li}) such that the vector of the
    nodes can be defined by \code{X = mu + li\%*\%E} where \code{E} is a
    normal random vector with expectation zero and variance matrix
    unity.
    \item  It is planned to add a fourth one under the name of /gbn/.
    
  }


  To relieve the memory effort, most names of the functions have been
  given a two (or more) components structure separated with a
  figure. This idea will be explained and exploited in a package to come
  named \pkg{documair}. The approximate meaning of the figures are:

  \itemize{

    \item 0 (similar to 'o') \code{rbmn0chain.01} to indicate an object
    example provided by \pkg{rbmn}.
    \item 1 (similar to an ~ and) ??? to link different objects or
    actions \code{train1car} for train and car.
    \item 2 (as usual but only one-to-one) \code{nbn2gema} means
    \"transforming a /nbn/ into a /gema/ objects\".
    \item 3 (remind the 'belong to' sign) \code{form3repeat} could be
    interpreted as "repeat action from the series of 'form' functions".
    \item 4 (associated to 'from') \code{adja4nbn} means "get the adjacency
    matrix from a /nbn/ object".
    \item 7 (upper bar of '7' similar to the hyphen) \code{arc7nb4nbn}
    means "get the arc-numbers from a /nbn/".
    \item 8 (similar to 'a') \code{generate8nbn} or \code{print8nbn}
    for \"generating or printing a /nbn/ object\".
  }

  A number of ancillary functions have not been exported to give a better
  access to the main function of /rbmn/. Nevertheless they are available
  in the \code{../rbmn/R/} directory, and with all their comments
  (equivalent to \code{Rd} files into \code{../rbmn/inst/original/}
  directory). Some of them are visible when defining the default
  arguments of some functions.

}

\section{Projected evolution of /mn/}{

  \itemize{
    \item Generalize the /mn/ object with a regression part like
    the output of function \code{condi4joint} when argument
    \code{pour} is not of length zero and argument \code{x2} is not
    null. With such a structure, every node of a /nbn/ could be
    described with a /mn/ comprising a unique variable... Also the two
    arguments of function \code{mn4joint1condi} would be just two /mn/
    objects... This is also the generalized /mn/ proposed in function
    \code{simulate8gmn} under the argument of \code{loi}... Of course
    almost all functions dealing with /nbn/ objects will be to rewrite!
    \item Introduce a new object \code{gbn} for Gaussian Bayesian
    network similar to the list provided by function \code{nbn2rr}.
  }

}

\section{TO DO list}{

  \itemize{
    \item Systemize the existence of \code{check8object} functions
    \item Introduce their systematic use conditionned with a
    \code{rbmn0check} variable.
    \item Follow the main checking of every functions
    \item Give (and use) class attributes to the main objects.
    \item Introduce the main objects in this short presentation.
    \item Make a true small example in this short presentation.
    \item Make the function \code{nbn4string7dag}.
    \item Add the computation made with /bnlearn/ in the example of
    \code{estimate8nbn}.
    \item Check the topological order within \code{nbn2nbn} depending on
    \code{rbmn0check} value.
    \item Make a super transformation function from an object associated
    to a Bayesian network to any other type, including itself.
    \item Correct the \code{ord} option in \code{order4chain}.
    \item Check the topological order in \code{rm8nd4adja}.
    \item Think about removing all \code{rmatrix} transformations to the
    benefit of the to-come \code{gbn} object.
    \item Introduce a check of non-negativity of \code{ma} into
    \code{cor4var}.
    \item Add examples to all functions without any.
  }
  
}

\author{

  Jean-Baptiste Denis\cr
  MIAj - Inra - Jouy-en-Josas \cr
  F-78532 Jouy-en-Josas

  Maintainer: Jean-Baptiste Denis \email{Jean-Baptiste.Denis@Jouy.Inra.Fr}

}



\references{

  (A technical report presenting the concepts used in \pkg{rbmn} is
  under redaction; it can be obtained as it is if asked.)

  Scutari M (2010). "Learning Bayesian Networks with the bnlearn R
      Package". Journal of Statistical Software, \strong{35}(3),
      1-22. URL http://www.jstatsoft.org/v35/i03/.

  Tian S, Scutari M & Denis J-B (2013, submitted to JSFdS). "Predicting
      with Crossed Linear Gaussian Bayesian Networks".
}
      
\examples{
library(rbmn)

## getting the data set
data(boco)
print(head(boco));

}
\keyword{ package }
