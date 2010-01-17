
# File NetworkEpiBayesSEIR.R

# FUNCTION itimestartvalues
# Internal function used by epibayesmcmc
# Determines start values for Infective times (I) if needed

itimestartvalues <- function(inf.list)
{	
	# Want to scroll through in order of increasing E times -- indexset gives this order
	indexset <- order(rank(inf.list[,3]))
	
	# Set initial I time
	inf.list[indexset[1],4] <- runif(1,inf.list[indexset[1],3],inf.list[indexset[2],3])
	currrectime <- nextrectime <- inf.list[indexset[1],5]
	nextrec <- indexset[1]
	
	# Cycle through the remaining individuals
	for (i in 2:length(inf.list[,4]))
	{		
		if (inf.list[indexset[i],3] > currrectime)
		{
			currrectime <- nextrectime
			inf.list[nextrec,4] <- runif(1,inf.list[nextrec,3],inf.list[indexset[i],3])
		}		
		inf.list[indexset[i],4] <- runif(1,inf.list[indexset[i],3],inf.list[indexset[i],5])
		if(inf.list[indexset[i],5] > nextrectime)
		{
			nextrectime <- inf.list[indexset[i],5]
			nextrec <- indexset[i]
		}
	}
	
	return(inf.list[,4])
}


# FUNCTION etimestartvalues
# Internal function used by epibayesmcmc
# Determines start values for Exposure times (E) if needed

etimestartvalues <- function(inf.list,initialoffset)
{	
	for (i in 1:length(inf.list[,3]))
	{
		currpar <- -999
		maxrand <- 0
		posspar <- 0
		for (j in 1:length(inf.list[,3]))
		{
			if (inf.list[i,4] > inf.list[j,4])	
			posspar <- runif(1,0,1)
			if (posspar > maxrand)
			{
				currpar <- j
				maxrand <- posspar
			}
		}
		if (currpar == -999) 
			inf.list[i,3] <- inf.list[i,4] - initialoffset
		else
			inf.list[i,3] <- runif(1,inf.list[currpar,4], min(inf.list[currpar,5],inf.list[i,4]))
	}	
	return(inf.list[,3])
}


# FUNCTION epibayesmcmc
# R wrapper for C function that performs Bayesian MCMC inference

epibayesmcmc <- function(inf.list, nsamp, thinning, bprior, tiprior, teprior, pprior = c(5,N), kiprior, 
	keprior, N=ninf, priordists="gamma", betapriordist=priordists, thetaipriordist=priordists, 
	thetaepriordist=priordists, kipriordist = priordists, kepriordist=priordists, extrathinning=FALSE,
	inferEtimes = FALSE, inferItimes = FALSE, parentprobmult = 1)
{		
	
	# Check overall input format
	
	if (!is.numeric(inf.list)) stop("Invalid input format: inf.list must by numeric")
			
	ninf <- length(inf.list[,5])
	
	if (N < ninf)	# number of nodes must be at least the number of infected
	{
		cat("Illegal N Value - setting N equal to number of infecteds!")
		N <- ninf
	}
	
	# Do some cleanup of parameters
	
	nsamp <- floor(nsamp)
	thinning <- floor(thinning)
	extrathinning <- floor(extrathinning)
	
	maxmove <- 11
	
	# Do some processing on infection/recovery times:
	
	# (1) Shift infection and recovery times so that the first recovery happens at time 0
	# 	 This is necessary since the prior on the initial infected is from -Inf to 0
	
	inf.list[,3:5] <- inf.list[,3:5] - min(inf.list[,5])
			
	# (2a) Set initial values of parameters
	# 	For beta, thetai, thetae, ki, ke, and p, we're starting these parameters at the mean of their respective prior distributions
	#		(except for the theta parameters -- we're starting them at their modes, since they may not have means)
	# 	Beta, ki, ke, thetai, and thetae have Gamma/IG (or uniform) priors, and p has a Beta prior
	
	if(betapriordist == "gamma") initbeta <- bprior[1]*bprior[2] else initbeta <- (bprior[1]+bprior[2])/2
	if(thetaipriordist == "gamma") initthetai <- tiprior[2]/(tiprior[1] + 1) else initthetai <- (tiprior[1]+tiprior[2])/2
	if(thetaepriordist == "gamma") initthetae <- teprior[2]/(teprior[1] + 1) else initthetae <- (teprior[1]+teprior[2])/2
	if(kipriordist == "gamma") initki <- kiprior[1]*kiprior[2] else initki <- (kiprior[1]+kiprior[2])/2
	if(kepriordist == "gamma") initke <- keprior[1]*keprior[2] else initke <- (keprior[1]+keprior[2])/2
	initp <- pprior[1]/(pprior[1]+pprior[2])

	# (2b) Set values for prior distributions
	#	0 indicates a uniform prior dist, and 1 indicates a gamma/IG prior dist (can expand this to more options later)
	#	p is assumed to have a beta dist in all cases
	
	bpriordistnum <- 1*(betapriordist == "gamma")
	tipriordistnum <- 1*(thetaipriordist == "gamma")
	tepriordistnum <- 1*(thetaepriordist == "gamma")
	kipriordistnum <- 1*(kipriordist == "gamma")
	kepriordistnum <- 1*(kepriordist == "gamma")
	
	# (3a) Check Removal time data
	
	if (sum(is.na(inf.list[,5])) > 0) stop("Invalid input data: NA not allowed for removal times.") 
	
	# (3b) Determine starting values for E/I times, if we're inferring these times (if not, make sure we have legitimate values)
	
	if (inferItimes) 
	{
		if (inferEtimes) 
			inf.list[,4] <- inf.list[,5] - rgamma(length(inf.list[,4]),initki,1/initthetai) 
		else 
		{
			if (sum(is.na(inf.list[,3])) > 0) stop("Invalid input data: NAs found in Exposure times -- run with InferEtimes = TRUE.") 
			if (sum(inf.list[,3] > inf.list[,5]) > 0) stop("Invalid input data: Exposure times cannot occur after Removal times.")
			inf.list[,4] <- itimestartvalues(inf.list) 
		}
	}
	else
		if (sum(is.na(inf.list[,4])) > 0) stop("Invalid input data: NAs found in Infective times -- run with InferItimes = TRUE.") 
		
	if (inferEtimes) 
	{
		if (sum(inf.list[,4] > inf.list[,5]) > 0) stop("Invalid input data: Infective times cannot occur after Removal times.")
		inf.list[,3] <- etimestartvalues(inf.list,initke*initthetae)	
	}
	else
		if (sum(is.na(inf.list[,3])) > 0) stop("Invalid input data: NAs found in Exposure times -- run with InferEtimes = TRUE.") 
	
	# (4) Set up storage for extra stuff returned (E/I times and transmission tree)
	#	Note: if we're not inferring the E/I times, then we just return NULL for these values
	#	There are generally too many times to return each iteration (we'll run out of space)
	#	So instead, if we want to return these, we use the variable "extrathinning" as the extra thinning interval
	
	if (extrathinning == 0)
	{
		storeexptimes <- NULL
		storeinftimes <- NULL
		storetranstree <- NULL
	} else
	{
		storeexptimes <- array(0,ninf*nsamp/(thinning*extrathinning))
		storeinftimes <- array(0,ninf*nsamp/(thinning*extrathinning))
		storetranstree <- array(0,ninf*nsamp/(thinning*extrathinning))
	}
		
	# (5) Do translation (renaming) of nodes and probable parent info
	inf.list[,1:2] <- floor(inf.list[,1:2])
	
	# (5a) Check for invalid parent references
	for (i in 1:length(inf.list[,2]))
		if (!is.na(inf.list[i,2]) & (!(inf.list[i,2] %in% inf.list[,1]) | (inf.list[i,2] == inf.list[i,1]))) 
			stop("Invalid input data: bad parent reference")
	
	# (5b) Fill in any missing IDs	
	nextid <-max(inf.list[,1],0,na.rm=TRUE)+1
	for (i in 1:length(inf.list[,1]))
		if (is.na(inf.list[i,1]))
		{
			inf.list[i,1] <- nextid
			nextid <- nextid +1
		}
		
	# (5c) Check for duplicate IDs	
	if(length(inf.list[,1]) > length(unique(inf.list[,1])))
		stop("Invalid input data: non-unique node IDs")
		
	# (5d) Rename Node IDs
	parentprior <- rep(0, times=length(inf.list[,2]))
	for (i in 1:length(inf.list[,2]))
		if(!is.na(inf.list[i,2]))
			parentprior[i] <- which(inf.list[,1] == inf.list[i,2])
		
	# Call C function to do actual MCMC routine
	
	output <- .C("epigraphmcmcc",as.double(inf.list[,3]),as.double(inf.list[,4]), as.double(inf.list[,5]),as.integer(nsamp),as.integer(thinning),as.double(bprior), 
		as.double(tiprior), as.double(teprior), as.double(pprior), as.double(kiprior), as.double(keprior), as.integer(ninf),as.integer(N), 
		as.double(initbeta),as.double(initthetai), as.double(initki), as.double(initthetae), as.double(initke), as.double(initp), as.integer(bpriordistnum), 
		as.integer(tipriordistnum), as.integer(tepriordistnum), as.integer(kipriordistnum), as.integer(kepriordistnum), accept=as.integer(array(0,maxmove-3)),
		propose=as.integer(array(0,maxmove-3)), betaout = as.double(array(0,nsamp/thinning)), thetaiout=as.double(array(0,nsamp/thinning)), 
		kiout=as.double(array(0,nsamp/thinning)), thetaeout=as.double(array(0,nsamp/thinning)), keout=as.double(array(0,nsamp/thinning)), 
		p=as.double(array(0,nsamp/thinning)),initexp=as.integer(array(0,nsamp/thinning)), initexptime=as.double(array(0,nsamp/thinning)), 
		exptimes = as.double(storeexptimes), inftimes = as.double(storeinftimes), transtree = as.integer(storetranstree), as.integer(extrathinning),
		as.integer(1*inferEtimes), as.integer(1*inferItimes), as.integer(parentprior), as.integer(parentprobmult),PACKAGE="epinet" 
		)
	
	# Post-processing	
	
	# Reverse translate nodes for initinf (back to original node IDs)
	initexpout <- inf.list[output$initexp,1]
	
	if(extrathinning == 0)
	{	
		expout <- NULL
		infout <- NULL
		transtreeout <- NULL
	} else
	{
		expout <- array(output$exptimes,dim=c(ninf,nsamp/(thinning*extrathinning)))
		infout <- array(output$inftimes,dim=c(ninf,nsamp/(thinning*extrathinning)))
		# Back translate transmission tree to their original node IDs
		output$transtree[output$transtree == -999] <- NA
		transtreeout <- array(inf.list[output$transtree,1],dim=c(ninf,nsamp/(thinning*extrathinning)))		
	}
	
	return(list(accept=output$accept,propose=output$propose,beta=output$betaout,thetai=output$thetaiout,thetae=output$thetaeout,ki=output$kiout,
		ke=output$keout, p=output$p,initexp=initexpout, initexptime=output$initexptime,exptimes=expout,inftimes=infout,rectimes=inf.list[,5],nodeid=inf.list[,1],
		transtree=transtreeout))
	
}
