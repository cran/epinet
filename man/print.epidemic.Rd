
\name{print.epidemic}

\Rdversion{2.10}

\alias{print.epidemic}

\title{Prints an epidemict object}

\description{Prints an object created by the
\code{\link{SEIR.simulator}} simulation routine.}

\usage{
	\method{print}{epidemic}(x, \dots)
	}

\arguments{
  \item{x}{an object of class epidemic, produced from the \code{\link{SEIR.simulator}} simulation function.}
  \item{\dots}{other arguments to be passed to the print routine.}
}

\details{ 
Prints the epidemic inference object, including the exposure, infectious, and recovery times
    of each node in the epidemic.
}

\value{ Strictly invoked for side effect.}

\seealso{\code{\link{SEIR.simulator}} for simulating an epidemic, \code{\link{summary.epidemic}}
for the summary method of an epidemic object, and \code{\link{plot.epidemic}} for plotting a
visual display of the epidemic.}

\references{
Groendyke, C. and Welch, D.  2018.  epinet: An R Package to Analyze Epidemics Spread across Contact Networks, \emph{Journal of Statistical Software}, \bold{83-11}.	}

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
print(exampleepidemic)
}

\keyword{graphs}
