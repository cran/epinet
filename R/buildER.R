
buildER <- function(N,p) 
  {
	  edgelist <- .Call("ERedgelist",as.integer(N),as.double(p),
  PACKAGE="epinet")
  return(matrix(data=edgelist,ncol=2,byrow=TRUE))
  }
  
