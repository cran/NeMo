\name{MixNet}
\alias{MixNet}

\title{MixNet graph generator}

\description{
  Generate a mixture of Erdos-Renyi random graph or \code{MixNet} graphs
  with \emph{q} classes. Erdos-Renyi graphs are generated from
  this function with one class.
}

\usage{
MixNet( N, Alpha, Pi, storage = "matrix", directed = FALSE,
        diag = FALSE )
}

\arguments{
  \item{N}{
    number of vertices.
    }
  \item{Alpha}{
    q-length vector of vertex proportions.
    }
  \item{Pi}{
    connectivity matrix with (\emph{q},\emph{q}) 
    dimensions. \code{Pi[i,j]} is the
    probability to connect from a vertex belonging to the j-class to a
    vertex belonging to the i-class. 
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

\value{ According to the \code{storage} parameter, returns an adjacency
  matrix or an edge list. For the adjacency matrix case,
  \code{A[i,j]=1} means that j is the source vertex and j
  is the destination vertex .  
}

\author{G.Grasseau}
\seealso{
 \code{\link{EDD}}
}

\examples{
#   Generate a Mixnet network of 100 vertices
#   with 3 classes
N <- 100

# Proportion vector of the 3 classes
Alpha <- c(0.3, 0.5, 0.2)

# Connectivity (probability) between the 3 classes
Pi <- matrix( c(
                 0.5 , 0.05, 0.1 ,
                 0.05, 0.4 , 0.01,
                 0.1 , 0.01, 0.3
               ), 3, 3 )

# Adjacency matrix
A <- MixNet (N, Alpha=Alpha, Pi=Pi )

# Edge list
A <- MixNet (N, Alpha=Alpha, Pi=Pi, storage="edges" )

}
\keyword{models}