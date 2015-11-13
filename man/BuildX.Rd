
\name{BuildX}

\alias{BuildX}

\title{Build a dyadic covariate matrix (X)}

\description{Build a dyadic covariate matrix (X) from a given nodal covariate matrix.}

\usage{BuildX(nodecov, unaryCol = NULL, unaryFunc = NULL, 
	binaryCol = NULL, binaryFunc = NULL, includeIntercept = TRUE) }

\arguments{
  \item{nodecov}{an N x k matrix where N is the number of nodes, column 1 is the node id and columns 2:k are covariate values for the node}
  \item{unaryCol}{a vector of column indices}
  \item{unaryFunc}{a vector of the same length as unaryCol of method names for comparing dyads.  Possible method names are "match" and "absdiff"}
  \item{binaryCol}{a list of 2 element vectors of column indices}
  \item{binaryFunc}{a vector of the same length as binaryCol of method names for comparing dyads.  Possible method names are "euclidean" and "manhattan"}
  \item{includeIntercept}{logical.  If TRUE, includes a column of all ones.  Defaults to TRUE}
}

\value{A dyadic covariate matrix with \eqn{{N \choose 2}}{(N choose 2)} rows, columns 1 and 2 are node ids, column 3 is all ones (if requested) and then one column for each given element of unaryCol and binaryCol.  

Assigns colnames depending on type of unaryFunc and binaryFunc and colnames of nodecov.
}

\seealso{\code{\link{SimulateDyadicLinearERGM}} for simulating a
    contact network based on a dyadic covairate matrix, and \code{\link{epinet}} for performing
    inference on the network and epidemic model parameters.
}


\author{  David Welch \email{david.welch@auckland.ac.nz}, Chris Groendyke \email{cgroendyke@gmail.com}
}

\examples{# make some nodal covariates
set.seed(3)
mycov = data.frame(id = 1:5, xpos = rnorm(5), ypos = rnorm(5), 
	house = c(1, 1, 2, 2, 2), gender = c(0, 0, 0, 1, 1))
# make matrix 
dyadCov = BuildX(mycov, unaryCol = c(4, 5), unaryFunc = c("match", "match"), 
	binaryCol = list(c(2, 3)), binaryFunc = "euclidean")
}

\keyword{graphs}
