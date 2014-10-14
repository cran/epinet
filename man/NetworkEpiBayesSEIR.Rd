
\name{epibayesmcmc}

\Rdversion{2.10}

\alias{epibayesmcmc}

\title{Uses epidemic data to perform Bayesian inference on a contact network}

\description{Performs Bayesian inference on parameters
for an SEIR epidemic model and a random graph
model, given recovery (and perhaps also exposure/infective) times
for each individual infected during the course of an epidemic.}

\usage{epibayesmcmc(epidata, dyadiccovmat, nsamp, thinning, bprior, tiprior, 
    teprior, etaprior, kiprior, keprior, etapropsd, 
    priordists = "gamma", betapriordist = priordists, 
    thetaipriordist = priordists,  thetaepriordist = priordists, 
    etapriordist = rep("normal", times=etapars), kipriordist = priordists, 
    kepriordist = priordists, extrathinning = FALSE,  inferEtimes = FALSE, 
    inferItimes = FALSE, parentprobmult = 1, verbose = TRUE)}

\arguments{
  \item{epidata}{input data consisting of exposure, infection, and recovery times.}
  \item{dyadiccovmat}{matrix of dyadic covariates.}
  \item{nsamp}{number of iterations to run MCMC algorithm.}
  \item{thinning}{thinning interval.}
  \item{bprior}{parameters for beta prior.}
  \item{tiprior}{parameters for thetai prior.}
  \item{teprior}{parameters for thetae prior.}
  \item{etaprior}{parameters for eta priors.}
  \item{kiprior}{parameters for ki prior.}
  \item{keprior}{parameters for ke prior.}
  \item{etapropsd}{standard deviation of proposal distributions for eta parameters.}
  \item{priordists}{can be \dQuote{uniform} or \dQuote{gamma}.}
  \item{betapriordist}{can be \dQuote{uniform} or \dQuote{gamma}.}
  \item{thetaipriordist}{can be \dQuote{uniform} or \dQuote{gamma}.}
  \item{thetaepriordist}{can be \dQuote{uniform} or \dQuote{gamma}.}
  \item{etapriordist}{prior distribution for the network parameters.}
  \item{kipriordist}{can be \dQuote{uniform} or \dQuote{gamma}.}
  \item{kepriordist}{can be \dQuote{uniform} or \dQuote{gamma}.}
  \item{extrathinning}{set to FALSE unless we want to return inferred values of the
  	exposure / infective times and the transmission tree, in which case it is an integer
  	specifying the extra thinning interval.}
  \item{inferEtimes}{boolean variable specifying whether the exposure times are considered 
 	 known (FALSE) or if they are unknown and need to be inferred (TRUE).}
  \item{inferItimes}{boolean variable specifying whether the infective times are considered 
 	 known (FALSE) or if they are unknown and need to be inferred (TRUE).}
  \item{parentprobmult}{multiplier for prior probability placed on suspected parent node.  Default
  	is a uniform prior assumption.}
  \item{verbose}{boolean variable specifying whether progress and information messages are
        displayed during the course of the MCMC routine.  Defaults to TRUE.}
}

\details{ Uses exposed, infective, and removal times from the infected nodes of an epidemic
in order to perform inference on the parameters of the network and epidemic models.   

epidata is an \eqn{N} row by 5 column array giving the identity, likely parent, and exposed, infective, 
and removal times for each of the \eqn{N} individuals in the population, as well as the values of any nodal 
covariates.  Column 1 gives the ID (an integer) of the node, and column 2 gives the identity of the probable 
parent of the node (if known).  Columns 3, 4, and 5 give the exposed, infective, and removal times.  Individuals who were not infected
during the course of the epidemic should have NA coded in columns 3, 4, and 5; the records for these
individuals should appear AFTER those corresponding to the individuals that were infected during the
epidemic.  Note that if the times are not internally consistent, an error message will be generated and no 
inference will be performed. It is necessary to include data for exposure and infective times, even 
if these values are not known (in this case, set the respective entries to NA).  

Any data rows corresponding to individuals not infected during the course of the epidemic, if present, 
must occur at the end of the array, after all rows for infected individuals.  These rows must have removal times of NA. 

dyadiccovmat is an \eqn{{N \choose 2}}{(N choose 2)} row by \eqn{(k+2)} matrix containing the dyadic covariates for the population, where \eqn{N} is the number of individuals in the population and \eqn{k} is the number of dyadic covariates used in the model.  The matrix contains one row for each dyad (pair of nodes).  Columns 1 and 2 give the ID of the two nodes comprising the dyad, and the remaining \eqn{k} columns give the covariate values.  Defaults to a matrix containing one covariate with a value of 1 for each dyad; this is the Erdos-Renyi model.

nsamp is the number of samples that will be produced for each of the model parameters. 

thinning is the thinning interval, e.g., to return every \eqn{10^{th}} sample, use thinning = 10.

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

The boolean variables inferEtimes and inferItimes indicate whether the exposure and infective times 
respectively are assumed to be unknown (need to be inferred by the algorithm) or are assumed to
be known (and are not updated in the algorithm).  If data exposure and / or infective times are not
known, set the respective columns of data to NA (any data passed to the function when will be ignored).

If exposure and / or infective times are being inferred and we wish to return the inferred values of these
times (along with the inferred transmission tree), set extrathinning equal to an integer specifying the
extra thinning interval for these values.  Because returning values for a large number of nodes can be
very space-intensive, an extra thinning interval can be given as a multiple of the thinning interval for
the other parameters.  For example, using thinning = 10 and extrathinning = 20 will return the values
of the inferred exposure and infective times and transmission tree every 200 iterations, and the values
of the other parameters every 10 iterations.  If these inferred values are not desired, set this variable to FALSE.

The default prior distribution for the parent of each node is uniform on all of the other nodes.  To specify a
non-uniform distribution, use column 2 of epidata and set parentpriormult to an integer multiplier greater than 1.

Uses an algorithm similar to that described in Groendyke et al. (2010) and Britton and O'Neill (2002).

}

\value{	
  \item{accept}{vector containing the number of times a proposed new value was accepted for the parameters 
      (P, eta, G, beta, thetai, ki, thetae, ke).}
  \item{propose}{vector containing the number of times a new value was proposed for the parameters 
      (P, eta, G, beta, thetai, ki, thetae, ke).}
  \item{llkd}{vector containing the log-likelihood at each iteration of the MCMC algorithm.}    
  \item{beta}{vector containing the sample for parameter beta.}
  \item{thetai}{vector containing the sample for parameter thetai.}
  \item{thetae}{vector containing the sample for parameter thetae.}
  \item{ki}{vector containing the sample for parameter ki.}
  \item{ke}{vector containing the sample for parameter ke.}
  \item{eta}{2-dimensional array containing the samples for the eta parameters.  The \eqn{i^{th}} column contains the 
  	sample for the \eqn{i^{th}} eta parameter.}
  \item{initexp}{vector containing the sample for parameter kappa (identity of initial exposed).
  	Will only vary when both the exposure and infection times are assumed unknown.}
  \item{initexptime}{vector containing the sample for parameter \eqn{E_{\kappa}}{E_kappa} (initial exposure time).
  	Will only vary when the exposure times are assumed unknown.}
  \item{exptimes}{if inferEtimes is TRUE, this is two-dimensional array containing the inferred exposure times
  	(exptimes[i, ] contains the sample of exposure times for node i). If inferEtimes is FALSE, this will be NULL.  }
  \item{inftimes}{if inferItimes is TRUE, this is two-dimensional array containing the inferred infective times
  	(inftimes[i, ] contains the sample of infective times for node i). If inferItimes is FALSE, this will be NULL.  }
  \item{rectimes}{vector containing the original recovery times.}
  \item{nodeid}{vector containing the node IDs for the individuals in the population.}
  \item{transtree}{A two-dimensional array containing the sample for inferred transmission tree.  transtree[i, ] contains the sample
  	of parent nodes for node i.  A parent node of -999 for i designates that i is the initial exposed node.}
}

\seealso{\code{\link{SEIR.simulator}} for simulating an SEIR epidemic over a contact network, and
    \code{\link{plotepitree}} for plotting the epidemic.} 

\references{
Groendyke, C., Welch, D. and Hunter, D.  2012. A Network-based Analysis of the 1861 
    Hagelloch Measles Data, \emph{Biometrics}, \bold{68-3}.

Groendyke, C., Welch, D. and Hunter, D.  2010. Bayesian inference for contact networks given
    epidemic data, \emph{Scandinavian Journal of Statistics}, \bold{38-3}.
	
Britton, T. and O'Neill, P.D. 2002. Bayesian inference for stochastic epidemics in populations 
    with random social structure, \emph{Scandinavian Journal of Statistics}, \bold{29-3}.
}

\author{ Chris Groendyke \email{cgroendyke@gmail.com},
    David Welch \email{david.welch@auckland.ac.nz} }

\examples{# Simulate an epidemic through a network of 30
set.seed(3)
N <- 30
# Build dyadic covariate matrix
# Have a single covariate for overall edge density; this is the Erdos-Renyi model
nodecov <- matrix(1:N,nrow = N)
dcm <- BuildDyadicCovMatrix(nodecov)
# Build network and simulate epidemic over network
examplenet <- BuildDyadicLinearERGM(N, dyadiccovmat=dcm, eta=-1.8)
exampleepidemic <- SEIR.simulator(examplenet, N = 30, 
    beta = 0.3, ki = 2, thetai = 5, latencydist="gamma")

# Run MCMC algorithm on this epidemic
# Note: Not enough data or iterations for any real inference
examplemcmc <- epibayesmcmc(exampleepidemic, nsamp = 5000, 
    thinning = 10, bprior = c(0, 1), tiprior = c(1, 3), teprior = c(1, 3), etaprior = c(0, 10),
    kiprior = c(2, 8), keprior = c(2, 8), etapropsd = 0.2, priordists = "uniform", 
    extrathinning = 10, inferEtimes = TRUE, inferItimes = TRUE)

}

\keyword{graphs}
