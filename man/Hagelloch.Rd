\name{Hagelloch}

\docType{data}

\alias{Hagelloch}

\title{Hagelloch measles data.}

\description{Epidemic data derived from a measles outbreak in the town of 
    Hagelloch, Germany in 1861.  188 individuals were infected over the 
    course of the epidemic.  Not all fields in the original data set have been 
    included in this package; the fields included here for each affected 
    individual are: NodeID (a numerical label); Putative Parent (node most 
    likely to have caused the infection, as recorded in the data); Exposure 
    time (time at which the individual contracted measles -- assigned NA 
    since these times are not given in the data); Infective time (time at which 
    the individual became infective -- these times were not directly given in 
    the data, but were assumed to be one day prior to the onset of prodromes); 
    and Removal time (time at which the individual ceased being infective) -- 
    these times were not directly given in the data, but were assumed to be 
    three days after the appearance of rash (or death, if sooner).  All times are 
    given in days.}

\usage{Hagelloch}

\format{Array consisting of 5 columns, in the format
expected by \code{\link{epibayesmcmc}}.  Column 1
gives the identity of the individual and column 2 gives
the identity of the putative parent node.  
Columns 3, 4, and 5 give the exposure, infective, 
and removal times for the individual.
}

\source{Thanks to Peter Neal for providing this data set.  This data was 
    originally collected by Pfeilsticker:
    
    Pfeilsticker, A. (1863). Beitrage zur Pathologie der Masern mit besonderer 
    Berucksichtigung der statistischen Verhaltnisse, M.D. thesis, Eberhard-Karls 
    Universitat, Tubingen.}

\references{Neal, P. and Roberts, G. (2004). Statistical inference and model selection for the 1861 
Hagelloch measles epidemic. \emph{Biostatistics} \bold{5} (2), 249.}

\keyword{datasets}

