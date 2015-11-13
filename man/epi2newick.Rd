
\name{epi2newick}

\Rdversion{2.10}

\alias{epi2newick}
\alias{epi2newickmcmc}

\title{Prints a transmission tree in Newick format.}

\description{Prints a simulated or inferred transmission tree in Newick format.}

\usage{epi2newick(epi)

epi2newickmcmc(mcmcoutput, index = dim(mcmcoutput$transtree)[2])}

\arguments{
  \item{epi}{a simulated epidemic, in the form of the
	output produced by \code{\link{SEIR.simulator}}.}
  \item{mcmcoutput}{output from \code{\link{epinet}}.}
  \item{index}{a number indicating which of the MCMC
	samples to plot. Defaults to the final sample in the chain.}
}

\details{ 
Converts the epinet epidemic format into a transmssion tree represented as a Newick string which is the standard tree format used in phylogenetics.  There are many packages available to analyse Newick format trees such as the ape package, IcyTree and FigTree. 
}

\value{A character string representing the epidemic transmission tree in Newick format.  Note that this string contains control characters that can be removed by using \code{\link{cat}}}

\seealso{\code{\link{epinet}} for generating posterior samples of the parameters, 
	\code{\link{print.epinet}} and \code{\link{summary.epinet}} for printing basic
	summary information about an epinet object, \code{\link{write.epinet}} for
	writing parameter and transmission tree posterior samples to file, and
    \code{\link{plotepitreemcmc}} for plotting the posterior samples of the transmission tree.} 

\references{
Rambaut A. 2014. FigTree v1.4. \url{tree.bio.ed.ac.uk/software/figtree/}.
Vaughan T. 2015. IcyTree \url{http://tgvaughan.github.io/icytree/icytree.html}.
}

\author{ David Welch \email{david.welch@auckland.ac.nz},
    Chris Groendyke \email{cgroendyke@gmail.com} }

\examples{
# Simulate an epidemic through a network of 30
set.seed(3)
N <- 30
# Build dyadic covariate matrix (X)
# Have a single covariate for overall edge density; this is the Erdos-Renyi model
nodecov <- matrix(1:N, nrow = N)
dcm <- BuildX(nodecov)
# Simulate network and then simulate epidemic over network
examplenet <- SimulateDyadicLinearERGM(N, dyadiccovmat = dcm, eta = -1.8)
exampleepidemic <- SEIR.simulator(examplenet, N = 30, 
    beta = 0.3, ki = 2, thetai = 5, latencydist="gamma")
cat(epi2newick(exampleepidemic))

\dontrun{
# Build covariates
set.seed(1)
N <- 50
mycov <- data.frame(id = 1:N, xpos = runif(N), ypos = runif(N))
dyadCov <- BuildX(mycov,binaryCol = list(c(2, 3)),binaryFunc = c("euclidean"))
# Build network
eta <- c(0, -7)
net <- SimulateDyadicLinearERGM(N = N,dyadiccovmat = dyadCov,eta = eta)
# Simulate epidemic
epi <- SEIR.simulator(M=net,N=N,beta=1,ki=3,thetai=7,ke=3,latencydist="gamma")
# Run MCMC routine on simulated epidemic
mcmcinput <- MCMCcontrol(nsamp = 1000000, thinning = 100, etapropsd = c(1, 1))
priors <- priorcontrol(bprior = c(0, 4), tiprior = c(1, 15), teprior = c(1, 15), 
	etaprior = c(0, 10, 0, 10), kiprior = c(1, 7), keprior = c(1, 7), priordists = "uniform")
out <- epinet(~ xpos.ypos.L2Dist, epidata = epi, dyadiccovmat = dyadCov,
	mcmcinput = mcmcinput, priors = priors)
cat(epi2newickmcmc(out))}
}

\keyword{graphs}
