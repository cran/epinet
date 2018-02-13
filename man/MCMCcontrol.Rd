
\name{MCMCcontrol}

\Rdversion{2.10}

\alias{MCMCcontrol}

\title{Set control parameters for epinet MCMC algorithm}

\description{Sets parameter values that control the MCMC
algorithm used by epinet to produce posterior samples.}

\usage{MCMCcontrol(nsamp, thinning, extrathinning = FALSE, burnin = 0, 
	seed = floor(runif(1, 0, 2^30)), etapropsd)}

\arguments{
  \item{nsamp}{number of iterations to run MCMC algorithm.}
  \item{thinning}{thinning interval.}
  \item{extrathinning}{set to FALSE unless we want to return inferred values of the
  	exposure / infective times and the transmission tree, in which case it is an integer
  	specifying the extra thinning interval.  Defaults to FALSE.}
  \item{burnin}{number of burn-in iterations to the run the MCMC algorithm.  Defaults to 0.}
  \item{seed}{seed for random number generation.  Defaults to a random value.}
  \item{etapropsd}{standard deviation of proposal distributions for eta parameters.}
}

\details{Auxiliary function that can be used to set parameter values that control the MCMC
algorithm used by epinet to produce posterior samples.  This function is only used in conjunction
with the \code{\link{epinet}} function.

nsamp is the number of samples that will be produced for each of the model parameters. 

thinning is the thinning interval, e.g., to return every \eqn{10^{th}} sample, use thinning = 10.

If exposure and / or infective times are being inferred and we wish to return the inferred values of these
times (along with the inferred transmission tree), set extrathinning equal to an integer specifying the
extra thinning interval for these values.  Because returning values for a large number of nodes can be
very space-intensive, an extra thinning interval can be given as a multiple of the thinning interval for
the other parameters.  For example, using thinning = 10 and extrathinning = 20 will return the values
of the inferred exposure and infective times and transmission tree every 200 iterations, and the values
of the other parameters every 10 iterations.  If these inferred values are not desired, set this variable to FALSE.

burnin controls the number of burn-in iterations to be run by the MCMC algorithm before samples begin
to become recorded.

etapropsd is a vector of length \eqn{k}, where \eqn{k} is the number of eta (network) parameters in the model,
including the intercept.  These are tuning parameters for the MCMC algorithm.
}

\value{A list with arguments as components.}

\seealso{ \code{\link{epinet}} for generating posterior samples of the parameters, and
    \code{\link{priorcontrol}} for specifying prior distributions and their hyperparameters.} 

\references{
Groendyke, C. and Welch, D.  2018.  epinet: An R Package to Analyze Epidemics Spread across Contact Networks, \emph{Journal of Statistical Software}, \bold{83-11}.	
}

\author{ Chris Groendyke \email{cgroendyke@gmail.com}}

\examples{# Simulate an epidemic through a network of 30
set.seed(3)
N <- 30
# Build dyadic covariate matrix (X)
# Have a single covariate for overall edge density; this is the Erdos-Renyi model
nodecov <- matrix(1:N, nrow = N)
dcm <- BuildX(nodecov)
# Simulate network and then simulate epidemic over network
examplenet <- SimulateDyadicLinearERGM(N, dyadiccovmat=dcm, eta=-1.8)
exampleepidemic <- SEIR.simulator(examplenet, N = 30, 
    beta = 0.3, ki = 2, thetai = 5, latencydist="gamma")
# Set inputs for MCMC algorithm
mcmcinput <- MCMCcontrol(nsamp = 5000, thinning = 10, etapropsd = 0.2) 
priorcontrol <- priorcontrol(bprior = c(0, 1), tiprior = c(1, 3), teprior = c(1, 3), 
    etaprior = c(0, 10), kiprior = c(2, 8), keprior = c(2, 8), priordists = "uniform")
# Run MCMC algorithm on this epidemic
# Note: Not enough data or iterations for any real inference
examplemcmc <- epinet( ~ 1, exampleepidemic, dcm, mcmcinput, priorcontrol)

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
}

}

\keyword{graphs}
