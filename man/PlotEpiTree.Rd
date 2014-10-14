

\name{PlotEpiTree}
\alias{plotepitree}
\alias{plotepitreemcmc}

\title{Plot the spread of an epidemic}

\description{Plot the spread of an epidemic over a contact network.}

\usage{plotepitree(epi, lwd = 1, leaf.labs = TRUE, leaf.cex = 0.75, 
    zero.at.start = FALSE, main = "Spread of Epidemic", xlab = "Time", 
    ylab= "", e.col = "black", i.col = "red", lty.transmission = 3, 
    marktransitions = TRUE, label.trans = "|", cex.trans = 0.5, ...)

plotepitreemcmc(mcmcoutput, index = lastiteration, lwd = 1, 
    leaf.labs = TRUE, leaf.cex = 0.75, zero.at.start = FALSE, 
    main = "Spread of Epidemic", xlab = "Time", ylab= "", 
    e.col = "black", i.col = "red", lty.transmission = 3, 
    marktransitions = TRUE, label.trans = "|", cex.trans = 0.5, ...)
}

\arguments{
		\item{epi}{a simulated epidemic, in the form of the
		    output produced by \code{\link{SEIR.simulator}}.}
		\item{mcmcoutput}{output from \code{\link{epibayesmcmc}}.}
		\item{index}{a number indicating which of the MCMC
		    samples to plot. Defaults to the final sample in the chain.}
		\item{lwd}{line width for the (horizontal) 
		    line segments representing the exposed and infective 
		    periods for each individual.  Also controls the line width 
		    for the (vertical) line segments showing the transmission
		    pathways.}
		\item{leaf.labs}{boolean variable controlling whether the leaf
		    labels (Node IDs) are displayed to the right of their
		    infective period.  Defaults to TRUE.}
		\item{leaf.cex}{Character expansion factor for the leaf
		    labels, if they are displayed.  Defaults to 0.75.}
		\item{zero.at.start}{boolean variable governing whether the
		    epidemic times are shifted so that the epidemic begins
		    at time zero.  Defaults to FALSE.}
		\item{main}{main title for plot.}
		\item{xlab}{label for x axis on plot.  Defaults to \dQuote{Time}.}
		\item{ylab}{label for y axis on plot.  Defaults to blank.}
		\item{e.col}{color to be used to plot the individuals' 
		    exposed periods on the plot.  Defaults to black.}
		\item{i.col}{color to be used to plot the individuals' 
		    infective periods on the plot.  Defaults to red.}
		\item{lty.transmission}{line type used to mark the (vertical)
		    infection pathway on the plot.  Defaults to 3 (dotted).}
		\item{marktransitions}{boolean variable indicating whether
		    tick marks should be placed at the times where the
		    individuals move from the exposed to the infective state.
		    Defaults to TRUE.}
		\item{label.trans}{character used to mark transition points, if
		    marktransitions is TRUE.  Defaults to \dQuote{|}.}
		\item{cex.trans}{magnification to be used to for transition labels,
		    if marktransitions is TRUE.  Defaults to 0.5.}
		\item{...}{other plotting arguments to be passed through to the
		    call to plot() that draws the axes and titles.}
}

\details{Takes as input either a simulated epidemic, or 
    the output from an MCMC output (which itself may 
    infer the epidemic times) and plots them, indicating the 
    path that the infection took during the epidemic (the 
    transmission tree) and the times that each node entered 
    the Exposed, Infective, and Removed states.  The default 
    plotting parameter values work well for epidemics up to
    about 50 - 60 infecteds and requires at least 2 infecteds.
    
    Only works for full data, i.e., the transmission tree must 
    be fully specified and all times for infected individuals
    must be known.  Remove any records containing missing 
    information prior to calling this function.}

\value{ returns no value.  Strictly invoked for the plotting
	side effect.
}

\seealso{\code{\link{SEIR.simulator}} for producing simulated
    epidemics and \code{\link{epibayesmcmc}} for performing
    inference on the network and epidemic model parameters.
}

\author{David Welch \email{david.welch@auckland.ac.nz}, 
    Chris Groendyke \email{cgroendyke@gmail.com}}

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

# Plot the simulated epidemic
plotepitree(exampleepidemic)

# Run MCMC algorithm on this epidemic
# Note: Not enough data or iterations to do inference
examplemcmc <- epibayesmcmc(exampleepidemic, nsamp = 5000, 
    thinning = 10, bprior = c(0, 1), tiprior = c(1, 3), teprior = c(1, 3), etaprior = c(0, 10),
    kiprior = c(2, 8), keprior = c(2, 8), etapropsd = 0.2, priordists = "uniform", 
    extrathinning = 10, inferEtimes = TRUE, inferItimes = TRUE, verbose = FALSE)

# Plot starting state of epidemic from chain 
plotepitreemcmc(examplemcmc, index = 1)

# Plot final state of epidemic from chain 
plotepitreemcmc(examplemcmc)
}

\keyword{graphs}
