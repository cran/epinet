
\name{buildER}
\alias{buildER}
\title{Build an Erdos-Renyi graph object}
\description{Builds an undirected Bernoulli (Erdos-Renyi) graph 
object and returns the result in the form of an edgelist matrix.}

\usage{buildER(N, p)}

\arguments{
  \item{N}{number of nodes in the graph.}
  \item{p}{probability of presence of an edge between two nodes.}
}

\details{Creates a random graph with N nodes. Generates edges
according to a Bernoulli random variable. Each potential edge
will exist with probability p. Uses an algorithm from Batagelj
and Brandes (2005). Partially based on C code from Butts (2008)}

\value{A 2 column matrix contains a list of
  all edges in the generated graph.}

\references{
	\itemize{
	\item{Batagelj, V. and Brandes, U. 2005. Efficient Generation of
		Large random Networks, \emph{Physical Review E}, \bold{71-3}}
	\item{Butts, C. 2008. network: A Package for Managing Relational
		Data in R, \emph{Journal of Statistical Software}, \bold{24-2}}
	\item{Erdos, P. and Renyi, A. 1959. On random graphs,
	\emph{Publicationes Mathematicae}, \bold{6}, 290-297.}  
	}
}

\author{Chris Groendyke \email{cxg928@psu.edu}}

\examples{# Create a Bernoulli network with 100 nodes

net <- buildER(100,0.01) }

\keyword{graphs}
