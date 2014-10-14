
plotepitree <- function(epi, lwd = 1, leaf.labs = TRUE, leaf.cex = 0.75, 
	zero.at.start = FALSE, main = "Spread of Epidemic", xlab = "Time", ylab= "",
	e.col = "black", i.col = "red", lty.transmission = 3, marktransitions = TRUE, 
	label.trans = "|", cex.trans = 0.5, ...){
        
        # find total number of infecteds -- need at least 2 infecteds in the epidemic
        ninf <- nrow(epi)
		if (ninf < 2) stop("Plot error: need at least two infecteds to plot the epidemic.")

	# remove susceptibles
		susc <- NULL
		for (i in 1:ninf) if(sum(is.na(epi[i,3:5]))==3) susc <- c(susc,i)
		if (length(susc)>0) epi <- epi[-susc,]	
        ninf <- nrow(epi)
		if (ninf < 2) stop("Plot error: need at least two infecteds to plot the epidemic.")
		
	# make sure we have full data
	if(sum(is.na(epi[,3:5])) > 0) stop("Plot error: need full data to plot the epidemic -- run again with all times known.")
	if(sum(is.na(epi[,1])) > 0) stop("Plot error: missing entries in first (Node ID) column.")
	if(sum(is.na(epi[,2])) > 1) stop("Plot error: too many NA values in second (Parent) column -- should only have one initial infected node.")
	
        if (zero.at.start) epi[,3:5] = epi[,3:5] - min(epi[,3:5])
 	
	# sort the epidemic by increasing order of exposure
	epi <- epi[order(epi[ ,3]), ]
	
	# count the cumulative number of children each infective has
	nchild = array(0,max(epi[ ,1]))
	for (i in ninf:2) nchild[epi[i,2]] = nchild[epi[i,2]] + nchild[epi[i,1]] + 1
	
	ypos = array(0,max(epi[,1]))
	ypos[epi[1,1]] = 1
	
	# position of the current is position of its parent + gap for as yet unencountered direct offspring of the parent
	for (i in 2:ninf) {  
		ypos[epi[i,1]] = ypos[epi[i,2]] + nchild[epi[i,2]] - nchild[epi[i,1]]
		nchild[epi[i,2]] = nchild[epi[i,2]] - nchild[epi[i,1]] - 1
	}
	
	# set up axes and titles
        plot(x = c(min(epi[ ,3]), max(epi[ ,5])+2), y = c(0.9,ninf ), pch=" ", 
		xlab = xlab, ylab = ylab, main = main, yaxt = "n", bty = "n", ...)      
        
        for (i in 1: ninf)        {
                # plot Exposed portion of infection
                lines(x=c(epi[i,3],epi[i,4]), y=c(ypos[epi[i,1]],ypos[epi[i,1]]), lwd = lwd, col = e.col)
                # plot Infective portion of infection
                lines(x=c(epi[i,4],epi[i,5]), y=c(ypos[epi[i,1]],ypos[epi[i,1]]), lwd = lwd, col = i.col)
                # plot connection between parent and child
                lines(x=c(epi[i,3],epi[i,3]), y=c(ypos[epi[i,2]],ypos[epi[i,1]]), lwd = lwd, lty = lty.transmission)
        }

        # plot labels
        if (leaf.labs) text(epi[,5],ypos[epi[,1]],labels = epi[,1],pos = 4, offset = 0.25, cex = leaf.cex)

	# Mark transition points
	if (marktransitions)
		for (i in 1:ninf)
			text(epi[i,4],ypos[epi[i,1]], labels=label.trans, cex = cex.trans)

	}

	
plotepitreemcmc <- function(mcmcoutput, index = lastiteration, lwd = 1,
	leaf.labs = TRUE, leaf.cex = 0.75, zero.at.start = FALSE, main = "Spread of Epidemic", 
	xlab = "Time", ylab= "", e.col = "black", i.col = "red", lty.transmission = 3, 
	marktransitions = TRUE, label.trans = "|", cex.trans = 0.5, ...) {

	lastiteration <- length(mcmcoutput$exptimes[1, ])
	if (lastiteration == 0) stop("Error: need inferred exposure and infectious times for this function")
	epi <- cbind(mcmcoutput$nodeid, mcmcoutput$transtree[ ,index], mcmcoutput$exptimes[ ,index], 
		mcmcoutput$inftimes[ ,index], mcmcoutput$rectimes)

	plotepitree(epi, lwd = lwd, leaf.labs = leaf.labs, leaf.cex = leaf.cex, 
		zero.at.start = zero.at.start, main = main, xlab = xlab, ylab= ylab, e.col = e.col, 
		i.col = i.col, lty.transmission = lty.transmission, marktransitions = marktransitions, 
		label.trans = label.trans, cex.trans = cex.trans, ...)
	
}
