
\name{epibayesmcmc}

\Rdversion{2.10}

\alias{epibayesmcmc}

\title{Uses epidemic data to perform Bayesian inference on a contact network}

\description{Performs Bayesian inference on parameters
for an SEIR epidemic model and a Bernoulli random graph
model, given recovery (and perhaps also exposure/infective) times
for each individual infected during the course of an epidemic.}

\usage{epibayesmcmc(inf.list, nsamp, thinning, bprior, tiprior, teprior, 
    pprior = c(5,N), kiprior, keprior, N, priordists="gamma", 
    betapriordist=priordists, thetaipriordist=priordists, 
    thetaepriordist=priordists, kipriordist = priordists, 
    kepriordist=priordists, extrathinning=FALSE, 
    inferEtimes = FALSE, inferItimes = FALSE, parentprobmult = 1)}

\arguments{
  \item{inf.list}{input data consisting of exposure, infection, and recovery times, such as produced by \code{\link{SEIR.simulator}}.}
  \item{nsamp}{number of iterations to run MCMC algorithm.}
  \item{thinning}{thinning interval.}
  \item{bprior}{list of parameters for beta prior.}
  \item{tiprior}{list of parameters for thetai prior.}
  \item{teprior}{list of parameters for thetae prior.}
  \item{pprior}{list of parameters for p prior.}
  \item{kiprior}{list of parameters for ki prior.}
  \item{keprior}{list of parameters for ke prior.}
  \item{N}{number of people in the population.  Assumed to be fixed and known.}
  \item{priordists}{can be \dQuote{uniform} or \dQuote{gamma}.}
  \item{betapriordist}{can be \dQuote{uniform} or \dQuote{gamma}.}
  \item{thetaipriordist}{can be \dQuote{uniform} or \dQuote{gamma}.}
  \item{thetaepriordist}{can be \dQuote{uniform} or \dQuote{gamma}.}
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
}

\details{ Uses exposed, infective, and removal times from the infected nodes of an epidemic
in order to perform inference on the parameters of the network and epidemic models.   

inf.list is an \eqn{m} row by 5 column array giving the identity, likely parent, and exposed, infective, 
and removal times for each of the \eqn{m} nodes that were infected during the course of the epidemic.  
Column 1 gives the ID (an integer) of the node, and column 2 gives the identity of the probable 
parent of the node (if known).  Columns 3, 4, and 5 give the exposed, infective, and removal times.  
Note that if these times are not internally consistent, an error message will be generated and no 
inference will be performed. It is  necessary to include data for exposure and infective times, even 
if these values are not known (in this case, set the respective columns of data to NA).

nsamp is the number of samples that will be produced for each of the model parameters. 

thinning is the thinning interval, e.g., to return every \eqn{10^{th}} sample, use thinning = 10.

The type of prior distribution (default is gamma / inverse gamma) can be specified for all parameters
using priordists or for each parameter individually.  The prior distribution for p is always the beta
distribution (set both hyper-parameters of this prior to 1 to get a uniform distribution).  For the other
parameters, either uniform or gamma / inverse gamma priors can be chosen.  (The two theta
parameters use inverse gamma prior distributions, while the other parameters use gamma priors.)

The parameters of the prior distributions are given as lists of (two) hyper-parameters.  If the uniform
prior is being used for a parameter, then the hyper-parameters are the lower and upper limits of
the distribution.  If the gamma distribution is being used with parameters \eqn{c} and \eqn{d}, then the prior mean
is \eqn{c \cdot d}{c * d} and the prior variance is \eqn{c \cdot d^2}{c * d^2}.  If the inverse gamma distribution is being used with parameters
\eqn{c} and \eqn{d}, then the prior mean is \eqn{\frac{d}{c-1}}{d/(c-1)} and the prior variance is 
\eqn{\frac{d^2}{(c-1)^2 \cdot (c-2)}}{d^2 / ( (c-1)^2 * (c-2) )}.

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
non-uniform distribution, use column 2 of inf.list and set parentpriormult to an integer multiplier greater than 1.

Uses an algorithm similar to that described in Groendyke et al. (2010) and Britton and O'Neill (2002).
}

\value{	
  \item{accept}{vector containing the number of times a proposed new value was accepted for the parameters 
      (P, p, G, beta, thetai, ki, thetae, ke).}
  \item{propose}{vector containing the number of times a new value was proposed for the parameters 
      (P, p, G, beta, thetai, ki, thetae, ke).}
  \item{beta}{vector containing the sample for parameter beta.}
  \item{thetai}{vector containing the sample for parameter thetai.}
  \item{thetae}{vector containing the sample for parameter thetae.}
  \item{ki}{vector containing the sample for parameter ki.}
  \item{ke}{vector containing the sample for parameter ke.}
  \item{p}{vector containing the sample for parameter p.}
  \item{initexp}{vector containing the sample for parameter kappa (identity of initial exposed).
  	Will only vary when both the exposure and infection times are assumed unknown.}
  \item{initexptime}{vector containing the sample for parameter \eqn{E_{\kappa}}{E_kappa} (initial exposure time).
  	Will only vary when the exposure times are assumed unknown.}
  \item{exptimes}{if inferEtimes is TRUE, this is two-dimensional array containing the inferred exposure times
  	(exptimes[i, ] contains the sample of exposure times for node i). If inferEtimes is FALSE, this will be NULL.  }
  \item{inftimes}{if inferItimes is TRUE, this is two-dimensional array containing the inferred infective times
  	(inftimes[i, ] contains the sample of infective times for node i). If inferItimes is FALSE, this will be NULL.  }
  \item{rectimes}{vector containing the original recovery times.}
  \item{transtree}{A two-dimensional array containing the sample for inferred transmission tree.  transtree[i, ] contains the sample
  	of parent nodes for node i.  A parent node of -999 for i designates that i is the initial exposed node.}
}

\seealso{\code{\link{SEIR.simulator}} for simulating an SEIR epidemic over an Erdos-Renyi network, \code{\link{Hagelloch}} 
    for an example epidemic data set, and \code{\link{PlotEpiTree}} for plotting the epidemic.} 

\references{
Groendyke, C., Welch, D. and Hunter, D. (2010) Bayesian inference for contact networks given
    epidemic data, Technical Report, Department of Statistics, Pennsylvania State University, \bold{10-02}
    \url{http://www.stat.psu.edu/reports/2010/TR10-02.pdf}
	
Britton, T. and O'Neill, P.D. 2002. Bayesian inference for stochastic epidemics in populations 
    with random social structure, \emph{Scandinavian Journal of Statistics}, \bold{29-3}
}

\author{ Chris Groendyke \email{cxg928@psu.edu}, 
    David Welch \email{jdw21@stat.psu.edu} }

\examples{# Simulate an epidemic through a network of 30
set.seed(1)
examplenet <- buildER(N = 30, p = 0.1)
exampleepidemic <- SEIR.simulator(examplenet, N = 30, 
    beta = 0.3, ki = 2, thetai = 5, latencydist = "gamma")

# Run MCMC algorithm on this epidemic
examplemcmc <- epibayesmcmc(exampleepidemic, nsamp = 100000, 
    thinning = 100, bprior = c(0, 1), tiprior = c(0, 5), teprior = c(0, 5), 
    kiprior = c(0, 10), keprior = c(0, 10), N = 30, priordists = "uniform")
}

\keyword{graphs}
