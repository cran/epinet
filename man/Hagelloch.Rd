\name{Hagelloch}

\docType{data}

\alias{HagellochTimes}
\alias{HagellochDyadCov}

\title{Hagelloch measles data.}

\description{Epidemic data derived from a measles outbreak in the town of 
    Hagelloch, Germany in 1861.  188 individuals were infected over the 
    course of the epidemic.  (One individual was removed from this dataset.)
    
    Consists of two files: HagellochTimes and HagellochDyadCov.  These two
    files contain the data necessary to analyze the Hagelloch measles data.
    
    HagellochTimes contains 187 rows (one for each individual included) and
    5 columns: NodeID (a numerical index ranging from 1 to 187); Putative
    Parent (the ID of the individual considered most likely to have been
    responsible for infecting this person, as determined by Osterle); Exposure 
    time, Infectious time, and Removal time (the time index in days at which 
    the individual entered the Exposed, Infectious, and Removed states,
    respectively).  Note that Exposure times are not known for this data set.
    See references below for details regarding the determination of 
    Infectious and Removal times.
    
    HagellochDyadCov is a matrix of dyadic covariates corresponding to
    the individuals in the Hagelloch data set.  Contains one row for each
    dyad (pair of individuals) in the population.  The first two columns are the
    Node IDs for the two individuals in the dyad.  Third column is a column
    of all 1 values, used as a baseline or intercept term.  Columns 4, 5, and 6
    are indicator variables for whether the two individuals in the dyad are in
    the same household, are both in classroom 1, or both in classroom 2,
    respectively.  Column 7 is the household distance between the two
    individuals in the dyad, measured in units of 2.5m, Columns 8 and 9 are indicator
    variables based on whether both individuals in the dyad are male or
    female, respectively.  Column 10 is the age difference (in years) between
    the two individuals in the dyad. 
    }

\usage{Hagelloch}

\format{See above.}

\source{Thanks to Peter Neal for providing this data set.  This data was 
    originally collected by Pfeilsticker:
    
    Pfeilsticker, A. (1863). Beitrage zur Pathologie der Masern mit besonderer 
    Berucksichtigung der statistischen Verhaltnisse, M.D. thesis, Eberhard-Karls 
    Universitat, Tubingen.
    
    and later modified by Osterle:
    
    Oesterle, H. (1992). Statistiche Reanalyse einer Masernepidemie 1861 in 
    Hagelloch, M.D. Thesis, Eberhard-Karls Universitat, Tubingen.

    }

\references{Neal, P. and Roberts, G. (2004). Statistical inference and model selection for the 1861 
Hagelloch measles epidemic. \emph{Biostatistics} \bold{5} (2), 249.}

\keyword{datasets}

