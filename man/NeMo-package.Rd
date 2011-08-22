\name{NeMo-package}
\alias{NeMo-package}
\alias{NeMo}
\docType{package}
\title{
  This package deals with discovering exceptional subgraphs or
  motifs in a network.
}
\description{
  The \pkg{NeMo} deals with  discovering exceptional subgraphs or
  motifs in a network. 
  The method implemented here to detect such exceptional motifs
  (or subgraphs) doesn't require any Monte-Carlo
  or sampling simulations.

  The main idea is to approximate the motif
  count distribution to compare with (p-value calculation) the observed motif
  count in the network of interest. As 
  shown in [JCB2008] a well-suited distribution is obtained for the compound
  Poisson distribution. Two parameters are required for this
  distribution and are easily derived from the distribution mean and
  the variance.

  Avoiding Monte-Carlo or sampling simulations, the mean
  and variance of the counts are analytically computed under random graph
  models. We propose three possible graph models: Erdos-Renyi (\code{ER}),
  Expected Degree Distribution (\code{EDD}) and MixNet model, known also
  as the Stochastic Block Model (\code{MixNet} for Erdos-Renyi Mixture for
  Network).

 The motif counting algorithm takes in account only \emph{non-induced}
 motifs, i.e motifs included in the subgraph extracted from the network
 (at the opposite, the \emph{induced} motifs are motifs which match exactly
 the subgraph extracted from the network).

 The main functions in this package are:
 \itemize{
   \item \code{getExceptMotif} computes the p-value of a motif to occur
   more frequently in the input network than those obtained in a random
   graph model.
   \item \code{countMotif} counts the motif occurrences in the input network.
 }

 The user can directly obtain the mean and variance of a motif counts
 for the three different random graph models and then compute the p-value  :
 \itemize{
   \item \code{getMotifMoment} returns the two first moments of the
   motif count distribution.
   \item \code{computePValue} returns the p-value of the observed motif
   counts assuming the motif count distribution is a geometric Poisson
   or Polya-Aeppli distribution.
 }
 
 Two more functions are provided to generate networks based on the
 different random graph models:
 
 \itemize{
   \item \code{EDD} generates \code{EDD} graphs,
   \item \code{MixNet} generates \code{MixNet} and \code{ER} graphs.
 }

 The node number of motifs is not limited to count the occurences of
 user specific motifs (user \emph{motif list}). Nevertheless,
 counting the occurences of \emph{all motifs} with k nodes is limited
 to 8 (8 for undirected graphs
 and 5 for directed graphs). For performance reasons, all these motifs
 with a given size are stored in files. 
 }

\section{Acknowledgement}{
This work has been supported by the French Agence Nationale de la
Recherche under grant NeMo ANR-08-BLAN-0304-01.
}

\author{
\itemize{
    \item Stéphane Robin \email{Stephane.Robin@agroparistech.fr},
    \item Etienne Birmelé \email{Etienne.Birmele@genopole.cnrs.fr},
    \item Michel Koskas \email{Michel.Koskas@agroparistech.fr},
    \item Gilles Grasseau \email{Gilles.Grasseau@genopole.cnrs.fr},
    \item Sophie Schbath \email{Sophie.Schbath@jouy.inra.fr}.
  }
}
\references{
[JCB2008] Picard F, Daudin JJ, Koskas M, Schbath S, Robin S.
\emph{Assessing the exceptionality of network motifs.}
Journal of Computational Biology 2008 Jan-Feb; Vol. 15 Number pp.1-20.
\url{http://www.ncbi.nlm.nih.gov/pubmed/18257674}
}
\keyword{ package }

