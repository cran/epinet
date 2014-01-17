######################################################################
#
# zzz.r
#
# .First.lib is run when the package is loaded with library(epinet)
#
######################################################################

.onAttach <- function(lib, pkg){
    packageStartupMessage("Welcome to epinet. Type help(package=\"epinet\") for help.\n")
}
