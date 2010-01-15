

# Function for plotting the output of an epidemic as produced by SEIR.simulator()

plotepi <- function(epi)
{
	
	plot(x=c(min(epi[,3]),max(epi[,5])),y=c(0,max(epi[,1]) ),pch=" ",xlab="Time",ylab="Infecteds",main="Spread of Epidemic",
		sub="\"+\" indicates time of transition from Exposed to Infective")

	ninf <- length(epi[,1])
	for (i in 1: ninf)
	{
		lines(x=c(epi[i,3],epi[i,5]),y=c(epi[i,1],epi[i,1]))
		arrows(x0=epi[i,3],y0=epi[i,2],x1=epi[i,3],y1=epi[i,1],length=min(max(3/ninf,0.05),0.3),col=i)
	}
	
	points(x=epi[,4],y=epi[,1],pch="+")

}


# Function for plotting the output of an epidemic whose E/I times have been inferred
# mcmcoutput is the output from  R function epibayesmcmc()
# num is the number of the sample to be plotted - by default, it plots the last sample

plotepimcmc <- function(mcmcoutput,num=iterations)
{
	iterations <- length(mcmcoutput$exptimes[1,])
	ninf <- length(mcmcoutput$exptimes[,1])
	epi <- cbind(mcmcoutput$nodeid,mcmcoutput$transtree[,num],mcmcoutput$exptimes[,num],mcmcoutput$inftimes[,num],mcmcoutput$rectimes)

	plot(x=c(min(epi[,3]),max(epi[,5])),y=c(0,max(epi[,1]+1)),pch=" ",xlab="Time",ylab="Infecteds",main="Spread of Epidemic",sub=
		"\"+\" indicates time of transition from Exposed to Infective")
	
	# Note: Need to play around with plotting parameters a bit more	
		
	for (i in 1:ninf)
	{
		lines(x=c(epi[i,3],epi[i,5]),y=c(epi[i,1],epi[i,1]))
		arrows(x0=epi[i,3],y0=epi[i,2],x1=epi[i,3],y1=epi[i,1],length=min(max(3/ninf,0.05),0.3),col=i)
	}
	
	points(x=epi[,4],y=epi[,1],pch="+")

}


