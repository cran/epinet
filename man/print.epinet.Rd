
\name{print.epinet}

\Rdversion{2.10}

\alias{print.epinet}

\title{Print basic information about an epinet object}

\description{Prints some general information about an object created by the
\code{\link{epinet}} inference routine.}

\usage{
	\method{print}{epinet}(x, \dots)
	}

\arguments{
  \item{x}{an object of class epinet, produced from the \code{\link{epinet}} inference function.}
  \item{\dots}{other arguments to be passed to the print routine.}
}

\details{ 
Prints some basic information about an epinet inference object, including the call, network
parameters in the model, and number of iterations of the MCMC algorithm.
}

\value{ Strictly invoked for side effect.}

\seealso{\code{\link{epinet}} for generating posterior samples of the parameters, and
    \code{\link{plot.epinet}} for plotting the posterior samples of the transmission tree.}

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
# Set inputs for MCMC algorithm
mcmcinput <- MCMCcontrol(nsamp = 5000, thinning = 10, etapropsd = 0.2) 
priorcontrol <- priorcontrol(bprior = c(0, 1), tiprior = c(1, 3), teprior = c(1, 3), 
    etaprior = c(0, 10), kiprior = c(2, 8), keprior = c(2, 8), priordists = "uniform")
# Run MCMC algorithm on this epidemic
# Note: Not enough data or iterations for any real inference
examplemcmc <- epinet( ~ 1, exampleepidemic, dcm, mcmcinput, priorcontrol)
print(examplemcmc)

\dontrun{
# Note: This may take a few minutes to run.
set.seed(1)
N <- 50
mycov <- data.frame(id = 1:N, xpos = runif(N), ypos = runif(N))
dyadCov <- BuildX(mycov, binaryCol = list(c(2, 3)),binaryFunc = c("euclidean"))
# Build network
eta <- c(0,-7)
net <- SimulateDyadicLinearERGM(N = N, dyadiccovmat = dyadCov, eta = eta)
# Simulate epidemic
epi <- SEIR.simulator(M = net, N = N, beta = 1, ki = 3, thetai = 7, ke = 3, latencydist = "gamma")
# Run MCMC routine on simulated epidemic
mcmcinput <- MCMCcontrol(nsamp = 1000000, thinning = 100, etapropsd = c(1, 1))
priors <- priorcontrol(bprior = c(0, 4), tiprior = c(1, 15), teprior = c(1, 15), 
	etaprior = c(0, 10, 0, 10), kiprior = c(1, 7), keprior = c(1, 7), priordists = "uniform")
out <- epinet(~ xpos.ypos.L2Dist, epidata = epi, dyadiccovmat = dyadCov,
	mcmcinput = mcmcinput, priors = priors)
print(out)
}

}

\keyword{graphs}
