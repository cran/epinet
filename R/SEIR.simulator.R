SEIR.simulator <- function(M, N, beta, ki, thetai, ke = ki, thetae = thetai, latencydist = "fixed", 
    latencyperiod = 0)
{
  
  # First do some input format checking

  # The line of code below converts a "network object" (as defined in the 
  # network package) to edgelist format.  This might be useful in some cases, 
  # but would require the package "epinet" to depend on the package 
  # "network", so it's not clear that it's worth it -- at least currently.
  
  #if(is.network(M)) M <- as.matrix.network(M,matrix.type="edgelist")

  # M should be an edgelist matrix
  if(!is.matrix(M))
  	stop("Input error: Network M must be an edgelist matrix.")
	
  # Make sure M has the correct dimensions
  if ( (length(dim(M)) != 2) || (dim(M)[2] != 2) )
	stop("Input error: Network M must an a 2-dimensional edgelist matrix.")
  
  # Begin actual simulator
	
  init <- sample(1:N,1)  # Inital infected individual is chosen at random

  # Keep a list of all upcoming transition and recovery times (t.time[i] = r.time = NA if i is susceptible)
  t.time = array(dim = N)
  r.time = array(dim = N)

  # Generate a transition time and recovery time for initial infected
  t.time[init] <- ifelse( latencydist=="fixed", latencyperiod, rgamma(1, ke, scale = thetae) )  
  r.time[init] = rgamma(1, ki, scale = thetai) + t.time[init]

  nextrec = init  	# Keep track of who, of current infecteds, is next to recover

  inf.list <- matrix(c(init, NA, 0, t.time[init], NA), nrow = 1)   # Keep track of initial infection

  time <- cm.time <- t.time[init]
    
  nexttrans = init # Temporary value
  t.time[init] <- Inf # Temporary value
  
  # Update storage with initial infected node
  s.list <- (1:N)[-init]
  e.list <- NULL
  i.list <- init

  inf <- list(i.list)		
  susc <- list(s.list)	
  expo<- list(e.list)		
   
  for( i in 2:(N*3) )	# Maximum number of iterations is 3*N (1 infection ,1 transition from  exposed to infective, and 1 recovery for every node)
  {
	    s.list<-array(susc[[i-1]])					# Which nodes are susceptible	
	    i.list<-array(inf[[i-1]])					# Which are infective
	    si.ex <- ((M[,1] %in% i.list) & (M[,2] %in% s.list)) | ((M[,2] %in% i.list) & (M[,1] %in% s.list))  # Extract indices of SI pairs
	    n.si<-sum(si.ex)							# Number of SI pairs
	
	    # Draw waiting times for the next removal, infection, and transition
	    dwt <- ifelse( length(inf[[i-1]]) > 0, r.time[nextrec] - cm.time, Inf ) 
	    bwt <- ifelse( n.si!=0, rexp(1,beta*n.si), Inf)
	    twt <- t.time[nexttrans] - cm.time
	    
	    ewt <- min(bwt, dwt, twt, na.rm = TRUE)
	    time <- c(time, ewt)			# Increment time
	    cm.time <- cm.time + ewt		# Increment cumulative time
	
	    if (ewt == bwt) test <- "Infect" else if (ewt == dwt) test <- "removal" else test <- "transition"
	    
	    if (test == "Infect")	# Event is an infection
	    {								
		      is.pairs <- which(si.ex == 1)	    # Choose indices where the infection will happen
	      
		      # For this infection, find the new infected node and its parent
		      smp.ind <- ifelse( length(is.pairs) == 1, is.pairs, sample(is.pairs,1) )  
		      parentindex <- which(M[smp.ind, ] %in% i.list)
		      new.inf <- M[smp.ind, 3-parentindex]
		      parent <- M[smp.ind, parentindex]
	
		      # Generate a transition time for new infected
		      lat <- ifelse( latencydist == "fixed", latencyperiod, rgamma(1, ke, scale = thetae) )
		      
		      t.time[new.inf] <- cm.time + lat
	
		      # Update lists of susceptible, exposed and infecteds
		      susc <- append(susc, list(susc[[i-1]][-which(susc[[i-1]] == new.inf)]))
		      expo <- append(expo, list(c(expo[[i-1]], new.inf))) 
		      inf <- append(inf, list(inf[[i-1]]))
	      
		      inf.list <- rbind(inf.list,c(new.inf, parent, cm.time, NA, NA))
		      
		      nexttrans <- which(t.time == min(t.time, na.rm = TRUE))
	
	    } else if (test == "removal")	# Event is a removal
	   {								
		      if(i==2)	# Only infected dies out
		      {
				  inf.list[1,5] <- cm.time
			      break
		      } 	
		
		      new.rec <- nextrec
		
		      # Update lists of susceptible, exposed and infecteds
		      susc <- append(susc, list(susc[[i-1]]))
		      expo <- append(expo, list(expo[[i-1]]))
		      inf <- append(inf, list(inf[[i-1]][-which(inf[[i-1]]==new.rec)]))      
		
		      inf.list[which(inf.list[,1]==new.rec), 5] <- cm.time
		      
		      # Update recovery times and next infected to recover
		      r.time[nextrec] <- NA
		      if(length(inf[[i]]) > 0) 
				nextrec <- which( r.time == min(r.time, na.rm = TRUE) ) 
		      else if (length(expo[[i]]) > 0)
		      {
				nextrec <- which( t.time == min(t.time, na.rm = TRUE) )
				r.time[nextrec] <- Inf
		      }
	      
	    } else 		# Event is a transition
	    { 	 
		      new.trans <- nexttrans
			 
		      # Update lists of susceptible, exposed and infecteds	 
		      susc <- append(susc,list(susc[[i-1]]))
		      expo <- append(expo,list(expo[[i-1]][-which(expo[[i-1]]==new.trans)]))
		      inf <- append(inf,list(c(inf[[i-1]],new.trans)))
		
		      inf.list[which(inf.list[,1]==new.trans),4] <- cm.time
		      
		      # Update transition times, recovery times, and assign a recovery time to the new transition
		      t.time[nexttrans] <- NA      
		      nexttrans <-  which(t.time == min(t.time,na.rm = TRUE))
		      r.time[new.trans] <- cm.time + rgamma(1,ki, scale = thetai)
		      if (r.time[new.trans] < r.time[nextrec]) nextrec <- new.trans
	    }
	    if(length(inf[[i]]) + length(expo[[i]]) ==0){break}	# No more infectious or exposed members, so epidemic ends
  }		
  
  inf.list[,3:5] <- inf.list[,3:5] - min(inf.list[,5])

  if (length(s.list) > 0)
  {
      for (i in 1:length(s.list))
          inf.list <- rbind(inf.list, c(s.list[i], NA, NA, NA, NA))
  }      

  colnames(inf.list) <- c("Node ID", "Parent", "Etime", "Itime", "Rtime")
  
  return(inf.list)
}
