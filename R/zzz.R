#
# copyright (c) 2009, Gilles Grasseau <gilles.grasseau@genopole.cnrs.fr>
# Last Modified 01/09/09
# Licensed under the GNU General Public License version 2 
#
######################################################################

.onLoad <- function(lib, pkg){
    setMotifsPath()
    ## ehelp <- help(package="paloma")$info[[1]]
    ## cat(paste(substring(ehelp[4],first=16),"\n",
    ##          "Version ",substring(ehelp[2],first=16),
    ##          " created on ",
    ##           substring(ehelp[3],first=16),".\n", sep=""))
    packageStartupMessage("")
    packageStartupMessage("   NeMo package")
    packageStartupMessage("   http://nemo.ssbgroup.fr")
    packageStartupMessage("")

    # Set the path of the motif data base 
}
