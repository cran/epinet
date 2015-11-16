
\name{write.epinet}

\Rdversion{2.10}

\alias{write.epinet}

\title{Writes posterior samples from an epinet object to an output file}

\description{Outputs posterior samples of an object created by the
\code{\link{epinet}} inference routine; creates two output files.}

\usage{write.epinet(out, filename)}

\arguments{
  \item{out}{an object of class epinet, produced from the \code{\link{epinet}} inference function.}
  \item{filename}{the name of the output file.}
}

\details{ 
Writes two output files corresponding to the output produced from the \code{\link{epinet}} 
inference function.  The first is a .log file, containing the posterior samples from the 
epidemic parameters in tab delimited form.  [This .log file can be read by Tracer, which calculates summary
statististics and diagnostics, and displays trace plots, histograms, etc.]  The second file
(which is only written if the transmission trees are returned from the inference routine),
is a .trees file, containing the inferred transmission trees, output in Newick format. 
}

\value{ Strictly invoked for side effect.}

\seealso{\code{\link{epinet}} for generating posterior samples of the parameters, 
	\code{\link{print.epinet}} and \code{\link{summary.epinet}} for printing basic
	summary information about an epinet object, \code{\link{epi2newickmcmc}} for
	printing an inferred transmission tree to the screen in Newick format, and
    \code{\link{plotepitreemcmc}} for plotting the posterior samples of the transmission tree.} 

\references{
Rambaut A., Suchard M., Xie D., Drummond A.J. 2014. Tracer v1.6. \url{http://beast.bio.ed.ac.uk/Tracer}.
}

\author{ David Welch \email{david.welch@auckland.ac.nz},
Chris Groendyke \email{cgroendyke@gmail.com} }

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
\dontrun{write.epinet(examplemcmc, "examplemcmc")}

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
write.epinet(out, "SampleInferenceOutput")
}

}

\keyword{graphs}
