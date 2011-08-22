\name{motif.list}
\alias{motif.list}

\docType{data}

\title{First motif lists} 

\usage{data(motif.list)}

\description{%
  Contains all the possible undirected motifs from 1-vertex motifs
  to 5-vertex motifs.
}

\format{%
  The first \code{motif.list} index points to the motif list with
  \emph{k} vertices while the second one gives the adjacency matrix of
  the motif. Thereby, \code{motif.list[[3]][[2]]} get the second
  motif adjacency matrix among the motifs with 3 vertices.
}

\examples{

# Get the motif lists
data(motif.list)

# All motifs with 4 vertices 
m4 <- motif.list[[4]]

#   Compute the p-value of all 4-vertex motifs
#   (by default EDD model)

# Build first the absolute path of one network example file
file.name <- system.file("networks", "Hpylo.txt", package="NeMo")

pv <- getExceptMotif( file.name,  motif.list[[4]] )
}

\keyword{datasets}
