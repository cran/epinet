
\name{SEIR.simulator}

\alias{SEIR.simulator}

\title{Simulate an epidemic on a contact network}

\description{Simulate the spread of an epidemic across an
    (undirected) contact network.}

\usage{SEIR.simulator(M, N, beta, ki, thetai, ke = ki, thetae = thetai, 
    latencydist = "fixed", latencyperiod = 0) }

\arguments{
  \item{M}{an undirected network, given in edgelist matrix 
      format}
  \item{N}{the number of nodes in the contact network.}
  \item{beta}{the transmission rate of the virus across an edge
      of the network.}
  \item{ki}{the shape parameter for the removal process for
      the epidemic.}
  \item{thetai}{the scale parameter for the removal process
      for the epidemic.}
  \item{ke}{the shape parameter for the removal process for
      the epidemic.}
  \item{thetae}{the scale parameter for the removal process
      for the epidemic.}
  \item{latencydist}{type of latency period; can be \dQuote{fixed} 
      or \dQuote{gamma}.}
  \item{latencyperiod}{length of latency period, if using latencydist 
      == \dQuote{fixed}.  Ignored if latencydist == \dQuote{gamma}.  
      Set to 0 to get an SIR epidemic.}
}

\details{Takes as input an undirected network, given in
edgelist matrix format, which is the same format returned by
\code{\link{SimulateDyadicLinearERGM}}. Randomly chooses an initial infective
individual. The infection spreads randomly across edges in the 
network according to exponential infective periods with mean \eqn{\frac{1}{beta}}{1 / beta}.
An infective individual remains in the exposed state for a either
a fixed period of time given by latencyperiod or a time described by a 
gamma RV with parameters ke and thetae (mean = \eqn{ke \cdot thetae}{ke * thetae}, 
var = \eqn{ke \cdot thetae^2}{ke * thetae^2}).  After this exposed period, an infected person 
moves to the Infected state, at which point they can infect susceptible individuals.  
The infective individuals are removed after an infective period whose length is governed by 
a gamma RV with parameters ki and thetai (mean = \eqn{ki \cdot thetai}{ki * thetai}, 
var = \eqn{ki \cdot thetai^2}{ki * thetai^2}).  Once an individual is removed, they cannot be 
re-infected and cannot infect others.}

\value{matrix consisting of one row for each individual in the 
  population.  Each row contains (in columns 1 - 5, respectively):
  the node infected, the infecting node, the time of infection,
  the time of transition from exposed to infective, and the 
  time of removal.  The times are shifted so that the first removal
  occurs at time 0.  The rows corresponding to the susceptible 
  members of the population (i.e., the individuals that were not 
  infected during the course of the epidemic) are placed after 
  those for the infected individuals.
}

\seealso{\code{\link{SimulateDyadicLinearERGM}} for simulating an Erdos-Renyi
    contact network, \code{\link{epinet}} for performing
    inference on the network and epidemic model parameters,
    and \code{\link{plot.epidemic}} and \code{\link{epi2newick}} for plotting functions.
}


\author{Chris Groendyke \email{cgroendyke@gmail.com},
  David Welch \email{david.welch@auckland.ac.nz}, 
}

\examples{# Simulate an epidemic through a network of 30
set.seed(3)
N <- 30
# Build dyadic covariate matrix (X)
# Have a single covariate for overall edge density; this is the Erdos-Renyi model
nodecov <- matrix(1:N,nrow = N)
dcm <- BuildX(nodecov)
# Simulate network and then simulate epidemic over network
examplenet <- SimulateDyadicLinearERGM(N, dyadiccovmat = dcm, eta = -1.8)
exampleepidemic <- SEIR.simulator(examplenet, N = 30, 
    beta = 0.3, ki = 2, thetai = 5, latencydist = "gamma")
}

\keyword{graphs}
