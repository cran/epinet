\name{ess}

\alias{ess}

\title{Calculate the Effective Sample Size}

\description{Calculate the Effective Sample Size for a marginal posterior sample obtained via MCMC}

\usage{ess(x, ignoreBurnin = FALSE, burninProportion = 0.1)
}

\arguments{
		\item{x}{a numeric vector of length N assumed to be samples from a Markov chain}
		\item{ignoreBurnin}{ logical indictating whether or not the first burninProportion of vector x should be ignored}
		\item{burninProportion}{if ignoreBurnin == TRUE, the first burninProportion*length(x) samples are removed from x before the ess is calculated}
}

\details{Calculates the effective sample size of x based on an estimate of the lag autocorrelation function.  Details of the method are in Section 11.5 of Bayesian Data Analysis, Third Edition, 2013, Andrew Gelman, John B. Carlin, Hal S. Stern, David B. Dunson, Aki Vehtari, Donald B. Rubin.}

\value{ Returns the estimated effective sample size for the last (1-burninProportion) samples in x.
}

\references{Gelman, A., Carlin, J.B., Stern, H.S., Dunson, D.B., Vehtari, A., Rubin, D.B., 2013 \emph{Bayesian Data Analysis}, Third Edition, (Section 11.5), Boca Raton, Florida: CRC Press.
}


\author{David Welch \email{david.welch@auckland.ac.nz}, 
    Chris Groendyke \email{cgroendyke@gmail.com}}

\examples{
	set.seed(8)
	x <- runif(1000)
	# expect ESS of close to 900 as samples are iid
	ess(x, ignoreBurnin = TRUE)
	# no burnin to ignore so ess is actually close to 1000
	ess(x, ignoreBurnin = FALSE)
	
	# ESS is a rough measure at best
	ess(1:1000,ignoreBurnin = FALSE)	
}
