
\name{priorcontrol}

\Rdversion{2.10}

\alias{priorcontrol}

\title{Set prior distributions and hyperparameters for epinet MCMC algorithm}

\description{Sets the prior distributions and corresponding hyperparameters
to be used in the epinet MCMC algorithm.}

\usage{priorcontrol(bprior, tiprior, teprior, etaprior, kiprior, keprior, 
	priordists = "gamma", betapriordist = priordists, thetaipriordist = priordists, 
	thetaepriordist = priordists, etapriordist = "normal", kipriordist = priordists, 
	kepriordist = priordists, parentprobmult = 1)}

\arguments{
  \item{bprior}{parameters for beta prior.}
  \item{tiprior}{parameters for thetai prior.}
  \item{teprior}{parameters for thetae prior.}
  \item{etaprior}{parameters for eta priors.}
  \item{kiprior}{parameters for ki prior.}
  \item{keprior}{parameters for ke prior.}
  \item{priordists}{can be \dQuote{uniform} or \dQuote{gamma}.}
  \item{betapriordist}{can be \dQuote{uniform} or \dQuote{gamma}.}
  \item{thetaipriordist}{can be \dQuote{uniform} or \dQuote{gamma}.}
  \item{thetaepriordist}{can be \dQuote{uniform} or \dQuote{gamma}.}
  \item{etapriordist}{prior distribution for the network parameters.}
  \item{kipriordist}{can be \dQuote{uniform} or \dQuote{gamma}.}
  \item{kepriordist}{can be \dQuote{uniform} or \dQuote{gamma}.}
  \item{parentprobmult}{multiplier for prior probability placed on suspected parent node.  Default
  	is a uniform prior assumption.}
}

\details{ Auxiliary function that can be used to set prior distributions and parameter 
	values that control the MCMC algorithm used by epinet to produce posterior 
	samples.  This function is only used in conjunction with the \code{\link{epinet}} 
	function.

The type of prior distribution (default is gamma / inverse gamma) can be specified for all epidemic parameters (i.e., all parameters except the eta network parameters)
using priordists or for each parameter individually.  Either uniform or gamma / inverse gamma priors can be chosen.  (The two theta
parameters use inverse gamma prior distributions, while the other epidemic parameters use gamma priors.)

The parameters of the epidemic parameter prior distributions are given as vectors of (two) hyper-parameters.  If the uniform
prior is being used for a parameter, then the hyper-parameters are the lower and upper limits of
the distribution.  If the gamma distribution is being used with parameters \eqn{c} and \eqn{d}, then the prior mean
is \eqn{c \cdot d}{c * d} and the prior variance is \eqn{c \cdot d^2}{c * d^2}.  If the inverse gamma distribution is being used with parameters
\eqn{c} and \eqn{d}, then the prior mean is \eqn{\frac{d}{c-1}}{d/(c-1)} and the prior variance is 
\eqn{\frac{d^2}{(c-1)^2 \cdot (c-2)}}{d^2 / ( (c-1)^2 * (c-2) )}.

For the network parameters (the eta parameters), the only prior assumption currently implemented is a set of independent normal distributions.

etaprior contains the hyper-parameters for the prior distributions of the eta parameters.  This is a vector of \eqn{2k} values,
giving the mean and standard deviation of each distribution (i.e., the first two entries are the mean and
standard deviation of the prior distribution for the first eta parameter, the next two entries are the mean and
standard deviation of the prior distribution for the second eta parameter, etc.)

The default prior distribution for the parent of each node is uniform on all of the other nodes.  To specify a
non-uniform distribution, use column 2 of epidata and set parentpriormult to an integer multiplier greater than 1.
}

\value{A list with arguments as components.}

\seealso{\code{\link{epinet}} for generating posterior samples of the parameters and
    \code{\link{MCMCcontrol}} for specifying control parameters for the MCMC algorithm.} 

\references{
Groendyke, C. and Welch, D.  2018.  epinet: An R Package to Analyze Epidemics Spread across Contact Networks, \emph{Journal of Statistical Software}, \bold{83-11}.	}

\author{ Chris Groendyke \email{cgroendyke@gmail.com} }

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
}

\keyword{graphs}
