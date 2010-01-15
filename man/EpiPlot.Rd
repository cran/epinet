
\name{EpiPlot}
\alias{EpiPlot}
\alias{plotepimcmc}
\alias{plotepi}
\title{Plot the spread of an epidemic}
\description{Plot the spread of an epidemic over a contact network.}

\usage{
plotepi(epi)
plotepimcmc(mcmcoutput,num=iterations)
}

\arguments{
		\item{epi}{a simulated epidemic, in the form of the
		output produced by \code{\link{SEIR.simulator}}.}
		\item{mcmcoutput}{output from \code{\link{epibayesmcmc}}.}
		\item{num}{a number indicating which of the MCMC
			samples to plot. Defaults to the final sample in the chain.}
}

\details{Takes as input either a simulated epidemic, or
	the output from an MCMC output (which itself infers
	the epidemic times) and plots them, indicating the
	path that the infection took during the epidemic (the
	transmission tree) and the times that each node entered
	the Exposed, Infectious, and Removed states.}

\value{ returns no value.  Strictly invoked for the plotting
	side effect.
}

\author{Chris Groendyke \email{cxg928@psu.edu}}

\examples{# Simulate an epidemic through a network of 30
set.seed(1)
examplenet <- buildER(30,0.1)
exampleepidemic <- SEIR.simulator(examplenet,30,0.3,2,5,latencydist="gamma")
# Plot the simulated epidemic
plotepi(exampleepidemic)
# Run MCMC algorithm on this epidemic
examplemcmc <- epibayesmcmc(exampleepidemic,nsamp=100000,
	thinning=1000, bprior=c(0,1), tiprior=c(0,5),teprior=c(0,5),kiprior=c(0,10),
	keprior=c(0,10), N=30, priordists="uniform",extrathinning=10,inferEtimes=TRUE, 
	inferItimes=TRUE)
# Plot starting state of epidemic from chain 
plotepimcmc(examplemcmc,1)
# Plot final state of epidemic from chain 
plotepimcmc(examplemcmc)
}

\keyword{graphs}
