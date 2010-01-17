######################################################################
#
# zzz.r
#
# .First.lib is run when the package is loaded with library(epinet)
#
######################################################################

.First.lib <- function(lib, pkg){
   library.dynam("epinet", pkg, lib)
    cat("Welcome to epinet. Type help(package=\"epinet\") for help.\n")
}
