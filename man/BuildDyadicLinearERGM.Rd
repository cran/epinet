
\name{BuildDyadicLinearERGM}

\Rdversion{2.10}

\alias{BuildDyadicLinearERGM}

\title{Builds an ERGM network using given covariate values}

\description{Simulates a random ERGM network using a given matrix of covariate
values and a corresponding vector of parameter values.}

\usage{BuildDyadicLinearERGM(N, dyadiccovmat, eta)}

\arguments{
  \item{N}{number of individuals in the population.}
  \item{dyadiccovmat}{matrix of dyadic covariates.}
  \item{eta}{vector of parameters.}
}

\details{ dyadiccovmat is an \eqn{{N \choose 2}}{(N choose 2)} by \eqn{(k+2)} matrix containing the dyadic covariates for the population, where \eqn{N} is the number of individuals in the population and \eqn{k} is the number of dyadic covariates used in the model.  The matrix contains one row for each dyad (pair of nodes).  Columns 1 and 2 give the ID of the two nodes comprising the dyad, and the remaining \eqn{k} columns give the covariate values.  \eqn{\eta}{eta} is the vector of parameters corresponding to the covariates.

For this class of dyadic independence network, the probability of an edge between individuals \eqn{i} and \eqn{j} is \eqn{p_{\{i,j\} }}{p_{i,j}}, where 
\deqn{\log \left( \frac{p_{\{i,j\} }}{1-p_{\{i,j\} }} \right) = \sum_{k} \eta_k X_{\{i,j\},k}}{log(p_{i,j} / (1 - p_{i,j})) = eta * X_{i,j} }

More information about this type of model can be found in Groendyke et al. (2012).

}

\value{a network in edgelist matrix format}

\seealso{\code{\link{SEIR.simulator}} for simulating an SEIR epidemic over a network. } 

\references{ Groendyke, C., Welch, D. and Hunter, D.  2012. A Network-based Analysis of the 1861 
    Hagelloch Measles Data, \emph{Biometrics}, \bold{68-3}.
}

\author{ David Welch \email{david.welch@auckland.ac.nz},
Chris Groendyke \email{cgroendyke@gmail.com}  }

\examples{
# Construct a network of 30 individuals
set.seed(3)
N <- 30
# Build dyadic covariate matrix
# Have a single covariate for overall edge density; this is the Erdos-Renyi model
nodecov <- matrix(1:N,nrow = N)
dcm <- BuildDyadicCovMatrix(nodecov)
# Build network and simulate epidemic over network
examplenet <- BuildDyadicLinearERGM(N, dyadiccovmat=dcm, eta=-1.8)
}

\keyword{graphs}
