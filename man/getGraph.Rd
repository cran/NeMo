\name{getGraph}
\alias{getGraph}

\title{Transform network data to a NeMoGraph object} 

\description{
This function is mainly used to handle a graph from an input file.
The returned value is a NeMoGraph object. Adjacency matrices and edge
lists can also be given in input. 
}
\usage{
getGraph( graph, directed = FALSE, verbose = FALSE)
}

\arguments{
  \item{graph}{
    describes the network, \code{graph} can be an adjacency matrix, an edge
    list \emph{ matrix(2, number.of.vertices) },
    a file name containing an edge/arc list (one edge/arc per line). 
  } 
  \item{directed}{
    \code{TRUE} if the network and motif are directed (otherwise
    \code{FALSE}).
  }   
  \item{verbose}{
    gives information messages about the reading steps.  
  }   
}   
\value{
  Returns an object discribing the network and containing the different fields:
  \item{directed}{ \code{TRUE} if the network is directed (otherwise
    \code{FALSE}).}
  \item{n.nodes}{number of nodes in the graph}
  \item{node.names}{node names}
  \item{n.edges}{ number of edges}
  \item{edges}{ edge list \emph{matrix(2, number.of.vertices)}. For
    directed network, \code{net$edges[1, ]} are the source vertices and
    \code{net$edges[2, ]} are the destination vertices}
}

\author{G. Grasseau}


\seealso{
\code{\link{countMotif}}, \code{\link{getExceptMotif}}
}
\examples{
#
#   To handle the network file
#

# Build absolute path of one network example file
file.name <- system.file("networks", "Hpylo.txt", package="NeMo")

network <- getGraph(file.name) 
counts <- countMotif( network, 4 )

# Use the network object
# ...
}

\keyword{file}
