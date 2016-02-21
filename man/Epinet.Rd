
\name{epinet}

\Rdversion{2.10}

\alias{epibayesmcmc}
\alias{epinet}

\title{Uses epidemic data to perform Bayesian inference on a contact network}

\description{Performs Bayesian inference on parameters
for an SEIR epidemic model and a random graph
model, given recovery (and perhaps also exposure/infective) times
for each individual infected during the course of an epidemic.}

\usage{epinet(formula, epidata, dyadiccovmat, mcmcinput = MCMCcontrol(), 
	priors = priorcontrol(), verbose = TRUE)}

\arguments{
  \item{formula}{an object of class \code{\link[stats]{formula}} giving a symbolic description of the model to be fit.}
  \item{epidata}{input data consisting of exposure, infection, and recovery times.}
  \item{dyadiccovmat}{matrix of dyadic covariates (X).  Can be constructed using \code{\link{BuildX}}.}
  \item{mcmcinput }{list of control options for MCMC algorithm.  Can be constructed using \code{\link{MCMCcontrol}}.}
  \item{priors}{list of prior distributions and parameters.  Can be constructed using \code{\link{priorcontrol}}.}
  \item{verbose}{boolean variable specifying whether progress and information messages are
        displayed during the course of the MCMC routine.  Defaults to TRUE.}
}

\details{ Uses exposed, infective, and removal times from the infected nodes of an epidemic
in order to perform inference on the parameters of the network and epidemic models.   

The formula will consist of variables (column names) found in the dyadiccovmat parameter.  By default, the model will include an intercept term.

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

dyadiccovmat is an \eqn{{N \choose 2}}{(N choose 2)} row by \eqn{(k+2)} column matrix containing the dyadic covariates for the population, where \eqn{N} is the number of individuals in the population and \eqn{k} is the number of dyadic covariates used in the model.  The matrix contains one row for each dyad (pair of nodes).  Columns 1 and 2 give the ID of the two nodes comprising the dyad, and the remaining \eqn{k} columns give the covariate values. 

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
  \item{exptimes}{if exposure times are inferred and corresponding posterior samples are returned, this is two-dimensional 
  	array containing the inferred exposure times (exptimes[\eqn{i}, ] contains the sample of exposure times for node \eqn{i}). 
	Otherwise, this will be NULL.  }
  \item{inftimes}{if infection times are inferred and corresponding posterior samples are returned, this is two-dimensional 
  	array containing the inferred infection times (inftimes[\eqn{i}, ] contains the sample of infection times for node \eqn{i}). 
	Otherwise, this will be NULL.  }
  \item{rectimes}{vector containing the original recovery times.}
  \item{nodeid}{vector containing the node IDs for the individuals in the population.}
  \item{transtree}{A two-dimensional array containing the sample for inferred transmission tree.  transtree[\eqn{i}, ] contains the sample
  	of parent nodes for node \eqn{i}.  A parent node of -999 for \eqn{i} designates that \eqn{i} is the initial exposed node.  If the transmission tree 
	is not inferred and returned, this will be NULL.}
  \item{call}{the matched call.}
  \item{formula}{the formula used in the inference routine.}
  \item{mcmcinfo}{input settings for the MCMC chain} 
}

\seealso{\code{\link{BuildX}} for building a dyadic covariate matrix, 
    \code{\link{MCMCcontrol}} for specifying control parameters for the MCMC algorithm,
    \code{\link{priorcontrol}} for specifying prior distributions and their hyperparameters, 
    \code{\link{epi2newick}} and \code{\link{write.epinet}} for writing the output of the algorithm
    to file, and \code{\link{plot.epinet}} for plotting the posterior samples of the transmission
    tree.} 

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
# Build dyadic covariate matrix (X)
# Have a single covariate for overall edge density; this is the Erdos-Renyi model
nodecov <- matrix(1:N,nrow = N)
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
