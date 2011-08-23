\name{getExceptMotif}
\alias{getExceptMotif}

\title{Compute the p-value of one or severals motif (subgraphs) in a
  network under a random graph model.}
\description{
  Main function of this package, \code{getExceptMotif} calculates the p-value
  that a motif occurs more frequently than expected in a given
  random model.
  
  First, the  mean and variance of the
  motifs counts are calculated analytically under a random graph
  model. Three reference graph models are possible: Erdos-Renyi (\code{ER}),
  Expected Degree Distribution (\code{EDD}) and the Erdos-Renyi Mixture for
  Network model
  (\code{MixNet} known also as the Stochastic Block Model).
  Then, the motif count distribution is approximated by a compound
  Poisson distribution (more precisely a geometric Poisson
  distribution, known also as Polya-Aeppli distribution) which depends
  on two parameters. These parameters
  are analytically computed from the mean and the variance of counts
  under the  chosen graph model \code{ER}, \code{EDD},
  \code{MixNet} (see \link{getMotifMoment}). 

  The p-value is obtained from the observed motif count in the input
  network (see \link{countMotif}) and the approximated count
  distribution, the geometric Poisson distribution.

  For for more details, please, see ref. [JCB2008].
}

\usage{
getExceptMotif( graph, motif = NULL, directed = FALSE, model = "EDD",
             verbose = FALSE, ...)
}

\arguments{
  \item{graph}{
    describes the network, \code{graph} can be an adjacency matrix, an edge
    list \emph{matrix(2, number.of.vertices)},
    a file name containing an edge/arc list (one edge/arc per line and
    indexes must be strictly positive) or a
    NeMoGraph object. 
  }
  \item{motif}{%
  three kinds of values are possible :
  \itemize{%
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
  \item{model}{ choose a random graphical model:
    \describe{%
      \item{\code{"ER"}}{ Erdos-Renyi graphical model. By default, the
	model parameter \code{Pi} is computed from the the input network.
	The \code{Pi} can be forced, see \ldots{} arguments.
      }
      \item{\code{"EDD"}}{ Expected Degree Distribution graphical model,
	By default, the model parameter \code{deg.dist} (degree
	distribution) is computed from the the input network.
	The \code{deg.dist} can be forced, see \ldots{} arguments.
      }
      \item{\code{"MixNet"}}{ Erdos-Renyi
	Mixture for Network. The model parameters \code{Alpha} and
	\code{Pi} must be provided. See \ldots{} arguments.}
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
	  \item{\code{deg.dist}}{degree distribution vector given in
	    input of the
	    Degree Distribution model. The \code{deg.dist[i]} value
	    corresponds to the degree \code{i-1}. For the directed case,
	    \code{deg.dist} is a matrix where
	    \code{deg.dist[ , 1] } are the incoming arc distribution and
	    \code{deg.dist[ , 2] } are the outgoing arc distribution.
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

\value{
  The value returned is a list with the following fields:
    \item{\code{model}}{
      graphical model used (\code{ER}, \code{EDD}
      \code{MixNet}) }
    \item{\code{motifs}}{
      list of motifs i.e. list of adjacency matrices}
    \item{\code{counts}}{
      motif counts in the
      given network (vector - one value per motif)}
    \item{\code{mean}}{
      theoretical mean given the random graph
      model (vector - one value per motif)}
    \item{\code{sd}}{
      theoretical standard deviation given the random
      graph model (vector - one value per motif)}
    \item{\code{a}}{
      Polya-Aeppli (geometric Poisson distribution)
      parameter  (vector - one value per motif)}
    \item{\code{lambda}}{
      Polya-Aeppli or geometric Poisson distribution
      parameter (vector - one value per motif)}
    
    \item{\code{P.value}}{
      p-value of the motif (vector - one value per motif)}
    \item{\code{param}}{
      parameters of the random graphical model
      (\code{Alpha}, \code{Pi}, \code{deg.dist})
    }
}
\author{
 S. Robin, E. Birmele, M. Koskas, G. Grasseau, S. Schbath
}

\seealso{
 \code{\link{getMotifMoment}}
}

\references{
[JCB2008] Picard F, Daudin JJ, Koskas M, Schbath S, Robin S.
\emph{Assessing the exceptionality of network motifs.}
Journal of Computational Biology 2008 Jan-Feb; Vol. 15 Number pp.1-20.
\url{http://www.ncbi.nlm.nih.gov/pubmed/18257674}
}

\examples{

# Build absolute path of one network example file
file.name <- system.file("networks", "Hpylo.txt", package="NeMo")

# Compute the p-value of all 4-vertex motifs
# (by default EDD model)
pv <- getExceptMotif( file.name, 4 )

#
#     Use Mixet Model
#

# Proportion vector of the 3 classes
Alpha <- c(0.3, 0.5, 0.2)

# Connectivity (probability) between the 3 classes
Pi <- matrix( c(
                 0.5 , 0.05, 0.1 ,
                 0.05, 0.4 , 0.01,
                 0.1 , 0.01, 0.3
               ), 3, 3 )

pv <- getExceptMotif( file.name, 4, model="MixNet", Alpha=Alpha, Pi=Pi )

#
#   To handle the network file
#
network <- getGraph(file.name) 
pv  <- getExceptMotif( network, 4, model="MixNet", Alpha=Alpha, Pi=Pi )
}

\keyword{stats}

