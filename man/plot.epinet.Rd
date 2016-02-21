

\name{plot.epinet}

\Rdversion{2.10}

\alias{plot.epinet}

\title{Plot the spread of an epidemic}

\description{Plot the spread of an epidemic over a contact network.}

\usage{\method{plot}{epinet}(x, index = dim(x$transtree)[2],
    lwd = 1, leaf.labs = TRUE, leaf.cex = 0.75, zero.at.start = FALSE, 
    main = "Transmission Tree", xlab = "Time", ylab= "", 
    e.col = "black", i.col = "red", lty.transmission = 3,
    marktransitions = TRUE, label.trans = "|", cex.trans = 0.5, ...)
}

\arguments{
		\item{x}{output from \code{\link{epinet}}.}
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

\details{Plots the output from the \code{link{epinet}} function
    indicating the path that the infection took during the
    epidemic (the transmission tree) and the times that each
    node entered the Exposed, Infective, and Removed states.
    The default plotting parameter values work well for
    epidemics up to about 50 - 60 infecteds.  For a larger
    number of infecteds, it is recommended to use
    \code{\link[grDevices]{pdf}} and adjust plotting dimensions.
}

\value{ returns no value.  Strictly invoked for the plotting
	side effect.
}

\seealso{\code{\link{epinet}} for performing
    inference on the network and epidemic model parameters.
    The functions \code{\link{epi2newick}} and 
    \code{\link{write.epinet}} offer other options for outputting
    inferred parameter samples and transmission trees.
}

\author{David Welch \email{david.welch@auckland.ac.nz}, 
    Chris Groendyke \email{cgroendyke@gmail.com}}

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
mcmcinput <- MCMCcontrol(nsamp = 5000, thinning = 10, extrathinning = 10, 
    etapropsd = 0.2) 
priorcontrol <- priorcontrol(bprior = c(0, 1), tiprior = c(1, 3), teprior = c(1, 3), 
    etaprior = c(0, 10), kiprior = c(2, 8), keprior = c(2, 8), priordists = "uniform")
# Run MCMC algorithm on this epidemic
# Delete Exposure and Infection times; will infer these in the MCMC algorithm
exampleepidemic[, 3:4] <- NA
examplemcmc <- epinet( ~ 1, exampleepidemic, dcm, mcmcinput, priorcontrol,
     verbose = FALSE)

# Plot starting state of epidemic from chain 
plot(examplemcmc, index = 1)

# Plot final state of epidemic from chain 
plot(examplemcmc)
}

\keyword{graphs}
