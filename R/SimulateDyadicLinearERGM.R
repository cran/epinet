

SimulateDyadicLinearERGM <- function(N, dyadiccovmat, eta)
{

	# Note: Only performing an error check on the number of rows and columns
	# Not checking order of dyads or for duplication
	# For now this is probably OK; need to coordinate with function that builds the dyadic cov matrix

	invlogit <- function(y) return(exp(y)/(1+exp(y)))
	
	etapars <- length(eta)
	numdyads <- N*(N-1)/2
	if(dim(dyadiccovmat)[1] != numdyads) stop("Invalid Input.")
	if(dim(dyadiccovmat)[2] != (etapars)+2) stop("Invalid Input.")
	
	heads <- NULL
	tails <- NULL
	
	counter <- 1
	for (i in 1:(N-1))
		for (j in (i+1):N)
		{				
			logodds <- eta * dyadiccovmat[counter,-(1:2)]
			if (runif(1,0,1) < invlogit(sum(logodds)))
			{
				heads <- c(heads,i)
				tails <- c(tails, j)
			}
			counter <- counter + 1
		}
	edgelistmat <- NULL
	if (!is.null(heads)) edgelistmat <- matrix(data=c(heads,tails),byrow=FALSE, nrow = length(heads), ncol=2)
	return(edgelistmat)	
}
