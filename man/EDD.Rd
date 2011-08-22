\name{EDD}
\alias{EDD}

\title{Expected Degree Distribution (EDD) graph generator}

\description{
  Generates graphs from a degree distribution. The probability to
  connect (edges/arcs) the vertex i to the vertex j is proportional
  to the product of the \emph{degree(i) . degree(j)}. See ref. [JCB2008]
  for more details. For \code{directed} networks, the total number of incoming
  arcs must be equal to the  number of outgoing arcs. If it isn't the
  case, an error message is returned.
}

\usage{
EDD( N, deg.dist, storage="matrix", directed = FALSE, diag=FALSE )
}
\arguments{
  \item{N}{
    number of vertices.
  }
  \item{deg.dist}{
    Degree Distribution. A vector for undirected graphs with
    \code{deg.dist[i]} corresponds to the degree \code{i-1}.
    For the directed case, \code{deg.dist} is a matrix where
    \code{deg.dist[ , 1]} are the incoming arc distribution and
    \code{deg.dist[ , 2]} are the outgoing arc distribution.
  }
  \item{storage}{
    storage mode: \code{matrix} means that an adjacency matrix is returned.
    \code{edges} means that an edge list is returned.
  }
  \item{directed}{
    \code{TRUE} if the network and motif are directed (otherwise
    \code{FALSE}).
  }
  \item{diag}{
    set to \code{TRUE} allows loops in the returned graph.
  }
}

\value{ Depending on \code{storage} parameter return an adjacency
  matrix or and an edge list. For the adjacency matrix case,
  \code{A[i,j] = 1} means that j is the source vertex and j
  is the destination vertex.  
}

\author{G.Grasseau}

\seealso{
 \code{\link{MixNet}}
}

\references{
[JCB2008] Picard F, Daudin JJ, Koskas M, Schbath S, Robin S.
\emph{Assessing the exceptionality of network motifs.}
Journal of Computational Biology 2008 Jan-Feb; Vol. 15 Number pp.1-20.
\url{http://www.ncbi.nlm.nih.gov/pubmed/18257674}
}

\examples{
#   Generate a EDD network of 100 vertices
N <- 100

# Degree distribution vector
deg.dist <- c(0.3, 0.4, 0.2, 0.1)

# Adjacency matrix
A <- EDD(N, deg.dist )

# Edge list
A <- EDD(N, deg.dist=deg.dist, storage="edges" )

# Degree distribution matrix for directed graphs
deg.dist.d <- cbind(
 c(0.3, 0.4, 0.2, 0.1),
 c(0.3, 0.4, 0.2, 0.1)
)

# Adjacency matrix
A <- EDD(N, deg.dist.d, directed=TRUE )
}  