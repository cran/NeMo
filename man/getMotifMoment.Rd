\name{getMotifMoment}
\alias{getMotifMoment}

\title{Computes the two first moments of the motif count distribution
  under a random graph model.}

\description{
  Computes the mean and the second
  moment of the count distribution, given the random graph model (\code{ER},
  \code{EDD} and \code{MixNet}). The variance of the counts, used in
  \link{getExceptMotif} is obviously derived from this two
  values (for more, see ref. [JCB2008]).  
}
\usage{
getMotifMoment( N, motif = NULL, directed = FALSE, model = "EDD",
             verbose = FALSE, ...)
}

\arguments{
  \item{N}{
    number of vertices in the network.
  }
  \item{motif}{%
    two kinds of values are possible :
    \itemize{%
      \item{an adjacency matrix describing the motif.}
      \item{a list of motifs i.e. a list of adjacency matrices.}
    }
  }
  \item{directed}{
    \code{TRUE} if the network and motif are directed (otherwise
    \code{FALSE}).
  }
  \item{model}{ chooses a random graphical model:
    \describe{%
      \item{\code{"ER"}}{ Erdos-Renyi graphical model. The
	model parameter \code{Pi} must be provided, see \ldots{} arguments.
      }
      \item{\code{"EDD"}}{ Expected Degree Distribution graphical model,
	The model parameter \code{deg.dist} (degree
	distribution) must be provided, see \ldots{} arguments.
      }
      \item{\code{"MixNet"}}{ Erdos-Renyi
	Mixture for Network. The model parameters \code{Alpha} and
	\code{Pi} must be provided, see \ldots{} arguments.}
    }
  }
  \item{verbose}{
    gives information messages about the computing steps. 
  }
  \item{\dots}{ graphical model parameters
    \itemize{
      \item{\code{ER} model:
	\describe{%
	  \item{\code{Pi}}{scalar value used with "ER" model. Connectivity
	    value in the Erdos-Renyi graphics model. Probability to connect
	    two vertices in the graph.}
	}
      }
      \item{\code{EDD} model:
	\describe{%
	  \item{\code{deg.dist}}{degree distribution (vector). The
	    \code{deg.dist[i]} value
	    corresponds to the degree \code{i-1}. For the directed case,
	    \code{deg.dist} is a matrix where
	    \code{deg.dist[ , 1] } describes the incoming arc distribution and
	    \code{deg.dist[ , 2] } describes the outgoing arc distribution.
	  }
	}
      }
      \item{\code{MixNet} model with \emph{q} classes:
	\describe{%
	  \item{\code{Alpha}}{proportion vector of the \emph{q} classes
	    given in input of the MixNet model.}
	  \item{\code{Pi}}{connectivity matrix given in input of the MixNet model.
	    \code{Pi[i,j]} is the probability
	    to connect from a vertices belonging to the j-class to a
	    vertex belonging
	    to the i-class. \code{Pi} is a matrix (q,q) where q is the
	    number of classes.}
	}
      }
    }
  }
}

\author{
 S. Robin, G. Grasseau
}

\references{
[JCB2008] Picard F, Daudin JJ, Koskas M, Schbath S, Robin S.
\emph{Assessing the exceptionality of network motifs.}
Journal of Computational Biology 2008 Jan-Feb; Vol. 15 Number pp.1-20.
\url{http://www.ncbi.nlm.nih.gov/pubmed/18257674}
}

\examples{
#
#     Mixet Model
#

# Proportion vector of the 3 classes
Alpha <- c(0.3, 0.5, 0.2)

# Connectivity (probability) between the 3 classes
Pi <- matrix( c(
                 0.5 , 0.05, 0.1 ,
                 0.05, 0.4 , 0.01,
                 0.1 , 0.01, 0.3
               ), 3, 3 )
motif <- matrix( c( 0, 1, 1, 1, 0, 1, 1, 1, 0), 3, 3)

moments <- getMotifMoment( 1000, motif, model="MixNet", Alpha=Alpha, Pi=Pi )

# Compute the 2 moments for all 4-vertex motifs.
data( motif.list)
moments <- getMotifMoment( 1000, motif.list[[4]], model="MixNet",
                            Alpha=Alpha, Pi=Pi )

}
	       