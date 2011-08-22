\name{countMotif}
\alias{countMotif}

\title{Count motifs or subgraphs in a network.}

\description{
\code{countMotif} counts the occurrence number of a motif (or a motif
list) in a network. Only the \emph{non-induced} motifs are counted
i.e. motifs included in subgraph extracted from the network.
For instance,
the 3-vertices subgraph with the \code{V} shape (2 edges) occurs 3 times
in a 3-vertices clique (triangle shape with 3 edges). At the opposite,
no occurrence are found in case of \emph{induced} subgraphs (no \code{V}
shape motif match exactly the 3-vertices clique). A fast backtracking
algorithm is implemented to count all the occurrences of the non-induced
subgraph.  
}
\usage{
countMotif( graph, motif, directed = FALSE, verbose = FALSE)
}
\arguments{
  \item{graph}{
    describes the network, \code{graph} can be an adjacency matrix, an edge
    list (\emph{matrix(2, number.of.edges))},
    a file name containing an egde/arc list (one edge/arc per line and
    indexes must be strictly positive), or
    a NeMoGraph object. 
  }
\item{motif}{ three kinds of values are possible :
  \itemize{
    \item{an adjacency matrix describing the motif.}
    \item{a list of motifs i.e. a list of adjacency matrices.}
    \item{a \emph{motif size} (numeric value) specifying to count all motifs
      having \emph{motif size} vertices. In this case, the motif size
      is limited to 8 for undirected motifs and 5 for directed motifs.} 
    }

  }
  \item{directed}{
    \code{TRUE} if the network and motif are directed (otherwise
    \code{FALSE}).
}
\item{verbose}{
  gives information messages about the computing steps.  
%%     ~~Describe \code{verbose} here~~
}
}
\value{
  Returns a list containing the different fields
  \item{motifs}{adjacency matrices list of the requested motifs.} 
  \item{counts}{occurrence numbers (vector) of the requested motifs.} 
}
\author{S. Robin, E. Birmelé, G. Grasseau, M. Koskas}

\section{Warning }{ Because the number of possible motifs increases
  exponentially, exploring all the motifs with \emph{motif size} greater
  than 8 is not tractable.
}

\seealso{
 \code{\link{getGraph}}
}
\examples{

#
# Count all 3-vertices subgraph from an adjacency matrix
#

# Define a "small" network.
network  <-  matrix( c( 0, 1, 1, 1, 0, 1, 1, 1, 0), 3, 3)

# Count all  3-vertices subgraph
counts <- countMotif( network, 3 )


#
#   From a file
#

# Build absolute path of one network example file
file.name <- system.file("networks", "Hpylo.txt", package="NeMo")

# Count all  4-vertices subgraph
counts <- countMotif( file.name, 4 )

#
#   To handle the network file
#
network <- getGraph(file.name) 
counts <- countMotif( network, 4 )
}
\keyword{models}
