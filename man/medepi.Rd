\name{medepi}
\docType{data}
\alias{medepi}
\title{Simulated epidemic data.}
\description{Data from a simulated SEIR epidemic
spread across a network of 30 individuals.  27 of
these individuals were ultimately infected.}
\usage{data(ExampleEpidemics)}
\format{Array consisting of 5 columns, in the format
expected by \code{\link{epibayesmcmc}}.  Column 1
gives the identity of the individual and column 2 gives
the identity of the parent node (individual that infected
the node).  Columns 3, 4, and 5 give the exposure,
infective, and removal times for the individual.}
%\source{}
%\references{}
\keyword{datasets}
