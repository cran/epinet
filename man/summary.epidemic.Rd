
\name{summary.epidemic}

\Rdversion{2.10}

\alias{summary.epidemic}

\title{Summarize simulated epidemic}

\description{Prints a summary of an epidemic simulated by the
\code{\link{SEIR.simulator}} simulation routine.}

\usage{
	\method{summary}{epidemic}(object, \dots)
	}

\arguments{
  \item{object}{an object of class epidemic, produced from the \code{\link{SEIR.simulator}} function.}
  \item{\dots}{other arguments to be passed to the summary routine.}
}

\details{ 
Prints a summary of the simulated epidemic, including the number of individuals infected over the
course of the epidemic, the number remaining susceptible throughout the epidemic, the total size
of the population, and length of the epidemic.
}

\value{Strictly invoked for side effect.}

\seealso{\code{\link{SEIR.simulator}} for simulating an epidemic, and
    \code{\link{plot.epidemic}} for plotting the simulated epidemic.}

\author{ Chris Groendyke \email{cgroendyke@gmail.com},
    David Welch \email{david.welch@auckland.ac.nz} }

\examples{# Simulate an epidemic through a network of 30
set.seed(3)
N <- 30
# Build dyadic covariate matrix (X)
# Have a single covariate for overall edge density; this is the Erdos-Renyi model
nodecov <- matrix(1:N, nrow = N)
dcm <- BuildX(nodecov)
# Simulate network and then simulate epidemic over network
examplenet <- SimulateDyadicLinearERGM(N, dyadiccovmat = dcm, eta = -1.8)
exampleepidemic <- SEIR.simulator(examplenet, N = 30, 
    beta = 0.3, ki = 2, thetai = 5, latencydist = "gamma")
summary(exampleepidemic)
}

\keyword{graphs}
