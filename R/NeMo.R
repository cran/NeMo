
getGraph <- function ( graph, directed=FALSE, verbose=FALSE ) {
# get a NeMo graph structure:
# + filename
# + edge list  
# + Adj matrix
# + NeMoGraph class
  
  # Argument "graph"
  #  transform to an edge list
  if( is.character( graph ) ) {
    G  <- readAdjacency( graph, directed=directed, diag=FALSE, verbose)
  } else if( is.matrix( graph ) ) {
    if (dim(graph)[1] ==  dim(graph)[2] ) {
      # Matrix case
      edges     <- AdjMat2Edges( graph, directed=directed)
      n.nodes   <- dim( graph ) [1]
      nodenames <- as.character( as.character(1:n.nodes) )
    } else if( dim(graph)[1] != 2) {
      # Edge List
      edges <- graph
      n.nodes <- max( edges )
      nodenames <- as.character( as.character(1:n.nodes) )
    }
    
    dimnames( edges ) [[1]] <- c("src", "dst")
    G <- list( directed=directed, n.nodes = n.nodes, node.names = nodenames,
               n.edges =dim( edges )[2], edges = edges )
  } else if ( is.NeMoGraph( graph ) ) {
    G <- graph
  } else {
    cat( "Error: bad \'graph\' argument \n" ) 
    stop()
  }
  class(G) <- "NeMoGraph" 
  invisible(G)
}

is.NeMoGraph <- function(x) { if (class(x)=="NeMoGraph") TRUE else FALSE }

countMotif <- function( graph, motif, directed=FALSE, verbose=FALSE ) {

  if( ! is.NeMoGraph( graph ) ) {
    graph <- getGraph( graph, directed=directed,  verbose=FALSE)
  }
  invisible( count.motifs( graph, motif, directed, verbose) )
  
}

count.motifs <- function( G, motif, directed=FALSE, verbose=FALSE ) {
# G : NeMo graph structure
  
  # Argument "motif"
  # and count motifs
  if (verbose) cat("\nCounting all motifs ...",  "\n")

  count  <- vector()
  if( ! is.list( motif ) & ! is.matrix(motif) ) {
    # Find all motif of a given size
    k.nodes <- motif
    if( (k.nodes > 0) &
       ( (!directed & (k.nodes <= 8 ))  | (directed & (k.nodes <= 5 )) ) ) {
      # Get motifs and the occurence number
      res  <- countAllNonExactMotif( G$edges, k.nodes,
                                    directed )
      nb.motifs <- length( res ) / 2 

      motif <- list()
      for( m in 1:nb.motifs ) {
        count[m]   <- res[[ nb.motifs + m ]]
        motif[[m]] <- matrix( res[[m]], k.nodes, k.nodes ) 
      }
      
    } else {
      cat( "Error: bad motif size\n" )
      stop()
    }
  } else {
    # Find motifs given in a list or matrix
    if ( is.matrix( motif ) ) {
       motif <- list( motif )
    }

    if ( is.list( motif ) ) {
      for( m in 1:length(motif) ) {
        count [m] <- countNonExactMotif( G$edges,
                                         AdjMat2Edges( motif[[m]] ),
                                         directed )
      }
    } else {
      cat( "Error: bad motif parameter\n" )
      stop()
    }
    
  }
  invisible( list( motifs=motif, counts=count) )
}

countNonExactMotif <- function( G, m, directed=FALSE ) {
# Call C function
# G : edge list of the network
# m : edge list of the motif

  # Arguments G
  # Pb for matrix (2,2)
  if( (dim( G )[1] != 2)  ) {
    cat( "Error: edge.list not a matrix [2, N] \n") 
    stop()
  }

  n <- max( G )
  k <- max( m )

  G.edge.nbr <- dim( G ) [2]
  m.edge.nbr <- dim( m ) [2]
  
  # C arrays start at 0
  G <- G - 1
  m <- m - 1

  # Set an end value to these lists
  G <- cbind( G, c(-1, -1) )
  m <- cbind( m, c(-1, -1) )
  
  directed.int <- directed*1
  return(
         .Call( "countNonExactMotif_C", G, n, m, k, directed.int,
                PACKAGE="NeMo")
         )
}

countAllNonExactMotif <- function( G, k, directed=FALSE ) {
# Call C function
# G : edge list of the network
# k : motif size

  # Pb for matrix (2,2)
  if( (dim( G )[1] != 2)  ) {
    cat( "Error: edge.list not a matrix [2, N] \n") 
    stop()
  }

  n <- max( G )

  G.edge.nbr <- dim( G ) [2]
  
  # C arrays start at 0
  G <- G - 1

  # Set an end value to these lists
  G <- cbind( G, c(-1, -1) )
  
  directed.int <- directed*1
  return(
         .Call( "countAllNonExactMotif_C", G, n, k, directed.int,
                PACKAGE="NeMo")
         )
}


computeMixNetMoment <- function( N, Alpha, Pi, motif, directed=FALSE ) {
# Compute E, and E^2 (2 first moments) of a motif in a MixNet model
# Alpha, Pi : MixNet model parameters
# motif : edge list

  testMixNetParameters( Alpha, Pi )

  Q <- length( Alpha )
  directed.int <-  directed*1

  ## Test if motif is an edge list 
  if( dim( motif )[1] != 2 ) {
    cat("Error: motif not an edge list \n"  )
    return(NULL)
  }
  
  k <- dim(motif)[2]
  motif[] <- motif[] - 1
  
  res <- .Call( "computeMixNetMoment_C", N, Q, Alpha, Pi,
                motif, k, directed.int, PACKAGE="NeMo")
}

computeEDDMoment <- function( N, Degree, motif, directed=FALSE ) {
# Compute E, and E^2 (2 first moments) of a motif in a EDD model
# Degree : MixNet model parameters
#   Degree  in  -> Degree[,1]
#   Degree  out -> Degree[,2]
  
  testEDDParameters( Degree, directed )
  n.degree <- length( Degree )
  directed.int <-  directed*1

  if( directed ) {
    n.degree <- n.degree / 2
  }


  ## Test if motif is an edge list 
  if( dim( motif )[1] != 2 ) {
    cat("Error: motif not an edge list \n"  )
    return(NULL)
  }

  k <- dim(motif)[2]
  motif[] <- motif[] - 1
  res <- .Call( "computeEDDMoment_C", N, as.vector( Degree ), n.degree,
                motif, k, directed.int, PACKAGE="NeMo")
}

computeERMoment <- function( N, Pi, motif, directed=FALSE ) {

  testERParameters( Pi )

  directed.int <-  directed*1

  ## Test if motif is an edge list 
  if( dim( motif )[1] != 2 ) {
    cat("Error: motif not an edge list \n"  )
    return(NULL)
  }
  k <- dim(motif)[2]
  motif[] <- motif[] - 1
  
  res <- .Call( "computeERMoment_C", N, Pi,
                motif, k, directed.int, PACKAGE="NeMo")

}

getMotifMoment <- function( N, motif=NULL, directed=FALSE,
                           model='EDD', verbose=FALSE, ... ) {
# graph :
# motif
# + motif size
# + adj. matrix
# + list of adj. matrix

# Model parameters
#  Alpha, Pi : MixNet 
#  Pi        (optionnal) : Erdos-Renyi (ER)
#  deg.dist  (optionnal) : Expected Degree Distribution (EDD)
  
  # Model parameters
  Alpha		<- sub.param("Alpha"		, NULL	, ...)
  Pi		<- sub.param("Pi"		, NULL	, ...)
  deg.dist	<- sub.param("deg.dist"		, NULL	, ...)

  if( model == "MixNet" ) {
    testMixNetParameters( Alpha, Pi )
  } else if( model == "EDD") {
    if (! is.null( deg.dist) ) {
      testEDDParameters( deg.dist, directed )
    }
  } else if( model == "ER") {
    if ( ! is.null( Pi ) ) {
      testERParameters ( Pi )
    }
  } else {
    cat( "Error: unknown model\n")
    stop()
  }
  if( model=='ER' ) {
    if( is.null( Pi ) ) {
      cat( "Error: for 'ER' model, Pi is required \n")
        stop()
    }
  }
  if( model=='EDD' ) {
    if( is.null( deg.dist ) ) {
      cat( "Error: for 'EDD' model, deg.dist is required \n")
      stop()
    }
  }
  
  if( verbose ) 
      cat("Computing moments ...",  "\n")

  # Transform in motif list
  if ( is.matrix( motif ) ) {
    motif <- list( motif )
  }
  
  #
  #  Compute moments
  #
  
  n.nodes <- N
  
  moments <- matrix( 0, length(motif), 2 )
  
  for( m in 1:length(motif) ) {
    name <- paste( "motif", m)
    
    if( verbose ) 
      cat("Motif", m,  "\n")
    
    if( model=='ER' ) {
      moments[m, ] <- computeERMoment( n.nodes, Pi,
                                    AdjMat2Edges( motif[[m]],
                                               directed=directed),
                                    directed=directed )
    
    } else if( model=='EDD' ) {
      moments[m, ] <- computeEDDMoment( n.nodes, deg.dist,
                                     AdjMat2Edges( motif[[m]],
                                                   directed=directed),
                                     directed=directed )
    } else if ( model=='MixNet') {
      moments[m, ] <- computeMixNetMoment( n.nodes, Alpha, Pi,
                                          AdjMat2Edges( motif[[m]],
                                                       directed=directed),
                                          directed=directed )
    }
  }
  invisible( moments )
}

getExceptMotif <- function( graph, motif=NULL, directed=FALSE,
                    model='EDD', verbose=FALSE, ... ) {
# graph :
# + filename
# + edge list  
# + Adj matrix  

# motif
# + motif size
# + adj. matrix
# + list of adj. matrix

# Model parameters
#  Alpha, Pi : MixNet 
#  Pi        (optionnal) : Erdos-Renyi (ER)
#  deg.dist  (optionnal) : Expected Degree Distribution (EDD)
  
  # Model parameters
  Alpha		<- sub.param("Alpha"		, NULL	, ...)
  Pi		<- sub.param("Pi"		, NULL	, ...)
  deg.dist	<- sub.param("deg.dist"		, NULL	, ...)

  if( model == "MixNet" ) {
    testMixNetParameters( Alpha, Pi )
  } else if( model == "EDD") {
    if (! is.null( deg.dist) ) {
      testEDDParameters( deg.dist, directed )
      if( is.vector( deg.dist ) ) {
        # Transform to a matrix[ ,2]
        deg.dist <- cbind( deg.dist, deg.dist)
      }
    }
  } else if( model == "ER") {
    if ( ! is.null( Pi ) ) {
      testERParameters ( Pi )
    }
  } else {
    cat( "Error: unknown model\n")
    stop()
  }
  
  # get a NeMo graph structure:
  if( is.NeMoGraph( graph ) ) {
    G <- graph
  } else {
    G <- getGraph( graph, directed, verbose )
  }

  if( is.null(deg.dist) & (model == "EDD") ) {
    deg.dist   <- CalculateDegreeDistribution( G$edges, directed=directed )
  }

  # Argument "motif"
  # and count motifs
  
  Counts <- count.motifs( G, motif, directed, verbose )
  
  if( verbose ) 
      cat("Computing moments ...",  "\n")

  #
  #  Compute moments
  #
  
  n.nodes <- G$n.nodes
  
  PA     <- vector()
  a      <- vector()
  lambda <- vector()
  mean.th   <- vector()
  Var.th <- vector()
  sd.th <- vector()
  param1 <- NULL
  param2 <- NULL
  for( m in 1:length(Counts$motifs) ) {
    name <- paste( "motif", m)
    
    if( verbose ) 
      cat("Motif", m,  "\n")

    if( model=='ER' ) {
      if( is.null( Pi ) ) {
        Pi <- 2 * dim(G$edges)[2] / ( n.nodes * (n.nodes -1))
        if(  directed ) {
          Pi <- Pi /2
        }
      }
      param1 <- Pi
      moments.th <- computeERMoment( n.nodes, Pi,
                                     AdjMat2Edges( Counts$motifs[[m]],
                                                   directed=directed),
                                    directed=directed )

    } else if( model=='EDD' ) {
      param1 <- deg.dist
      if( ! directed ) {
        moments.th <- computeEDDMoment( n.nodes, deg.dist[ ,1],
                                       AdjMat2Edges( Counts$motifs[[m]],
                                                    directed=directed),
                                       directed=directed )
      } else {
        moments.th <- computeEDDMoment( n.nodes, deg.dist,
                                       AdjMat2Edges( Counts$motifs[[m]],
                                                    directed=directed),
                                       directed=directed )
      }
    } else if ( model=='MixNet') {
      param1 <- Alpha
      param2 <- Pi
      moments.th <- computeMixNetMoment( n.nodes, Alpha, Pi,
                                          AdjMat2Edges( Counts$motifs[[m]],
                                                       directed=directed),
                                          directed=directed )
    }
    
    mean.th  [m] <- moments.th[1]
    Var.th[m] <- moments.th[2] -  moments.th[1]^2
    if ( abs( Var.th[m] ) < 1.e-6 )
      Var.th[m] <- 0.0
    sd.th[m] <- sqrt( Var.th[m] )
    a     [m] <- (Var.th[m] - mean.th[m]) /
                   (Var.th[m] + mean.th[m])
    lambda[m] <- (1 - a[m])  * mean.th[m]
    
    if( verbose ) {
      cat("  mean.th, sd.th, count  :", mean.th[m], sd.th[m], Counts$counts[m],"\n")
      cat("  a, lambda, count  :", a[m], lambda[m],  Counts$counts[m],"\n")
    }
    
    if (verbose) cat("  Computing Polya-Aeppli ...",  "\n")

    # PA    [m] <- PolyaAeppli.uppertail(lambda[m], a[m], Counts$counts[m] )
    # PA    [m] <- pgeopoi( Counts$counts[m] , lambda[m], (1.0 - a[m]), lower.tail=FALSE)

    PA[m] <- computePValue( lambda[m], a[m], Counts$counts[m] )
    
    if (verbose) cat("  PA :", PA[m], "\n" )
  }

  if( verbose ) {
    cat(" Motif         Nobs               mean.th           sd.th      P-value   \n")
    cat(" ------------------------------------------------------------------------\n")
    for( m in 1:length(Counts$motifs) ) {
      str <- sprintf( " %5d %20d  %15.2f  %15.2f  %.2E \n",
                     m, Counts$counts[m], mean.th[m], sd.th[m] , PA[m] )
      
      cat( str )
    }
  }
  invisible( list(model=model, motifs=Counts$motifs, counts=Counts$counts,
                  mean=mean.th, sd= sd.th, a=a, lambda=lambda, P.value = PA,
                  param=list( param1,param2 ) )) 
}

# ----------------------------------------------------------- */
# ...Fonction qui calcule la queue de distribution à gauche...*/
# ...d'une loi de Polya-Aeppli de parametres (lambda, a) : ...*/
# ...P(PA<=n) avec n un entier                             ...*/
# ...D'apres Nuel (2006)                                   ...*/
# ----------------------------------------------------------- */

#ATTENTION : pb avec cette fonction en R :
#Error in if ((S.suiv > 0) && (S.suiv < INFINITY)) A.suiv <- A.cour else { : 
#missing value where TRUE/FALSE needed



PolyaAeppli.lowertail<-function(lambda, a, n) 
{
  INFINITY<-1e300
  PRECISION<-1e-15


  z<- (1.0-a) * lambda / a ; 

  # ........ Calcul de la queue inférieure P(PA <= n) : ................. */

  
  if (n==0)
    res<-exp(-lambda) 
  else if (n==1)
    res<-exp(-lambda)*(1.0+a*z)
  else
    {
      L.prec <- -lambda 
      L.cour <- L.prec + log(a*z) 
      A.cour <- L.prec 
      S.cour <- 1.0 + a*z 
      for (i in 2:n)
	{
	  L.suiv <- L.cour + log( a*(2.0*i-2.0+z)/i + a*a*(2.0-i)*exp(L.prec-L.cour)/i ) 
	  S.suiv <- S.cour + exp(L.suiv-A.cour) 
	  if ( (S.suiv>0) && (S.suiv<INFINITY))
	    A.suiv <- A.cour 
	  else
	    {
	      A.suiv <- A.cour + log(S.cour) 
	      S.suiv <- 1.0 + exp(L.suiv-A.suiv) 
	    }
	  
	  L.prec <- L.cour 
	  L.cour <- L.suiv 
	  S.cour <- S.suiv 
	  A.cour <- A.suiv 
	  
	} 
      res <- S.cour * exp(A.cour) 
    }
  res

}

# ----------------------------------------------------------- */
# ...Fonction qui calcule la queue de distribution à droite...*/
# ...d'une loi de Polya-Aeppli de parametres (lambda, a) : ...*/
# ...P(PA>=n) avec n un entier                             ...*/
# ...D'apres Nuel (2006)                                   ...*/
# ----------------------------------------------------------- */
PolyaAeppli.uppertail<-function(lambda, a, n) 
{
  INFINITY<-1e300
  PRECISION<-1e-15


  converged<-FALSE
  z <- (1.0-a) * lambda / a  
  # ........ Calcul de la queue supérieure P(PA >= n) : ................. */  

  # GG ??? 3 lines
  if( (a) < 0 ) {
    cat("Warning PA: log( negative ) \n" )
    res <- 0
  } else {

  if (n==0)
    res <- 1.0 
  else if (n==1)
    res <- 1.0 - exp(-lambda)
  else
    {
      L.prec <- -lambda   # ln P(PA=0) */
      L.cour <- L.prec + log(a*z) # ln P(PA=1) */
      for (i in 2:n)  #récurrence pour calculer ln P(PA=i) */
	{
	  L.suiv <- L.cour + log( a*(2.0*i-2.0+z)/i + a*a*(2.0-i)*exp(L.prec-L.cour)/i ) 
	  L.prec <- L.cour 
	  L.cour <- L.suiv      
	}
      A.cour <- L.suiv   # ln P(PA=n) */
      S.cour <- 1.0 
      # i=n+1 */
      while (converged==FALSE)
	{
	  L.suiv <- L.cour + log( a*(2.0*i-2.0+z)/i + a*a*(2.0-i)*exp(L.prec-L.cour)/i ) 

          # GG mofif
	  # if (L.suiv < log(PRECISION) + A.cour + log(S.cour))
          # cat( "L.suiv =", L.suiv, " A.cour =", A.cour, " log(S.cour) =" , log(S.cour), "\n")
	  if ( L.suiv < ( log(PRECISION) + A.cour + log(S.cour) ) )
	    {
	      converged <- TRUE
	    } 
	  
	  S.suiv <- S.cour + exp(L.suiv-A.cour) 
	  if ( (S.suiv>0) && (S.suiv<INFINITY))
	    A.suiv <- A.cour 
	  else
	    {
	      A.suiv <- A.cour + log(S.cour) 
	      S.suiv <- 1.0 + exp(L.suiv-A.suiv) 
	    } 
	  L.prec <- L.cour 
	  L.cour <- L.suiv      
	  S.cour <- S.suiv 
	  A.cour <- A.suiv   
	  i<-i+1
          # cat("S.cour = ",  S.cour, ", exp(A.cour) = ", exp(A.cour) ,"pv = ", S.cour * exp(A.cour), "\n")
	}
      
      res <- S.cour * exp(A.cour) 
    }
  }
  res

}


PolyaAeppli<-function(lambda,a,n)
{
  z <- (1.0-a) * lambda / a
  if (n==0)
    res <- exp(-lambda)
  else if (n==1)
    res <- a * z * exp(-lambda)
  else
    {
      p.prec<-exp(-lambda)
      p.cour<-a * z * exp(-lambda)
          
      for (i in 2:n)
        {
          p.suiv <- (2*i-2+z)*a*p.cour/i + (2-i)*a*a*p.prec/i
          p.prec<-p.cour
          p.cour<-p.suiv
       }
      res<-p.cour
    }
  res
}

#
# P-value computation
#
computePValue <- function( lambda, a, counts ) {
  if( (a >= 1) | (a <= 0) ) {
    cat("Error: bad Polya-Aeppli parameter, 0 < a < 1 \n"  )
    return(NA)
  }
  if( (lambda <= 0) ) {
    cat("Error: bad  Polya-Aeppli parameter, lambda > 0 \n"  )
    return(NA)
  }
  if( (counts < 0) ) {
    cat("Error: bad  Polya-Aeppli parameter, counts > 0 \n"  )
    return(NA)
  }

#  if( counts < ( mean.th +  sd.th ) ) {
#    pa <-  1.0 - PolyaAeppli.lowertail( lambda, a, counts-1 )
#  } else {
#    pa <-  1.0 - PolyaAeppli.uppertail( lambda, a, counts-1 )
#  }

  if( counts == 0 ) {
    pa <- 0.0 
  } else {
    pa <- PolyaAeppli.lowertail( lambda, a, counts-1 )
  }
  return( 1.0  - pa )
}


CalculateDegreeDistribution <- function( edge.list,  directed=FALSE ) {

  if ( ! is.matrix( edge.list ) ) {
    cat( "ERROR : edge.list not a matrix [2, N] \n") 
    stop()
  } else if ( dim( edge.list )[1] != 2) {
    cat( "ERROR : edge.list not a matrix [2, N] \n") 
    stop()
  }
  
  Adj <- Edges2AdjMat( edge.list, symmetrize=!directed )
  deg.out <- colSums( Adj )
  deg.in  <- rowSums( Adj )

  n.nodes <- length( deg.out  )
  max.deg <- max( deg.in , deg.out) + 1
  
  dist.deg.in  <- rep( 0, max.deg )
  dist.deg.out <- rep( 0, max.deg )
  for (i in 1:n.nodes ) {
    dist.deg.in [ deg.in [i] + 1 ] = dist.deg.in [ deg.in [i] + 1 ] + 1
    dist.deg.out[ deg.out[i] + 1 ] = dist.deg.out[ deg.out[i] + 1 ] + 1
  }
  dist.deg <- cbind(  dist.deg.in, dist.deg.out )
  return( dist.deg )
}

# Polya-Aeppli Gregory
dgeopoi=function(x,lambda,theta) {

res=rep(0,x+1);

z=lambda*theta/(1-theta);

res[1]=exp(-lambda);

res[2]=exp(-lambda)*(1-theta)*z;
exp

# GG Add test
if( x > 2 ) {
for (i in 3:(x+1) ){
    n=i-1;
    res[i]=(2*n-2+z)/n*(1-theta)*res[i-1]+(2-n)/n*(1-theta)^2*res[i-2];
}
}
res;
}


# lower.tail=TRUE return P(N<=x)
# lower.tail=FALSE return P(N>x)
pgeopoi=function(x,lambda,theta,lower.tail=TRUE) {
res=0;
x=round(x);
z=lambda*theta/(1-theta);
if (lower.tail==TRUE) {
   # compute P(N<=x)
   aux1=exp(-lambda);
   if (x==0) {
      res=aux1;
   } else {
     aux2=exp(-lambda)*(1-theta)*z;
     res=aux1+aux2;
     if (x>=2)
     for (n in 2:x) {
	 aux=(2*n-2+z)/n*(1-theta)*aux2+(2-n)/n*(1-theta)^2*aux1;
	 res=res+aux;
	 aux1=aux2;
	 aux2=aux;
     }
   }
} else {
  x=x+1;
  # compute P(N>=x)
  if (x==0) {
     res=1.0;
  }
  if (x==1) {
     res=ppois(x-1,lambda,lower.tail);
  }
  if (x>=2) {
#     y=exp(-lambda)*lambda*theta*(1-theta)^(x-1);
#     aux1=1/theta;
#     res=aux1*y;
#     aux2=((x-2)*theta+1)/(x-1)/theta^2;
#     y=y*z*(x-1)/2;
#     res=res+aux2*y;
#     if (x>=3)  
#     for (m in 3:x) {
#	 #print(c(m,aux1,aux2,y));
#	 aux=(x-m+2)*(x-m+1-(x-m+2-m)*(1-theta))*aux2+(x-m+1)*(2-m)*(1-theta)*aux1;
#	 aux=aux/theta/(x-m+1)/(x-m+2);
#	 y=y*z*(x-m+1)/(m-1)/m;
#	 res=res+aux*y;
#	 aux1=aux2;
#	 aux2=aux;
#     }
     # switch to log scale
     y=-lambda+log(lambda*theta*(1-theta)^(x-1));
     aux1=1/theta;
     res=aux1*exp(y);
     aux2=((x-2)*theta+1)/(x-1)/theta^2;
     y=y+log(z*(x-1)/2);
     res=res+aux2*exp(y);
     aux1=log(aux1);
     aux2=log(aux2);
     if (x>=3)  
     for (m in 3:x) {
	 #print(c(m,aux1,aux2,y));
	 aux=(x-m+2)*(x-m+1-(x-m+2-m)*(1-theta))+(x-m+1)*(2-m)*(1-theta)*exp(aux1-aux2);
	 aux=aux/theta/(x-m+1)/(x-m+2);
	 y=y+log(z*(x-m+1)/(m-1)/m);
	 res=res+exp(aux2+y)*aux;
	 aux1=aux2;
	 aux2=aux2+log(aux);
     }
  res=res+ppois(x,lambda,lower.tail);
  } 
}
res;
}


#postscript("geopoi_pdf.eps");
#
#x=50; mean=10;
#main="lambda/theta=10.0 theta=1.0,0.9,...,0.1 (from top to bottom in x=10)";
#xlabel="x";
#ylabel="P(N=x)";
#res=dpois(0:x,lambda=mean);
#plot(0:x,res,main=main,xlab=xlabel,ylab=ylabel,t='l',col=0);
#theta=0.9;
#res=dgeopoi(x,lambda=mean*theta, theta=theta);
#lines(0:x,res,t='l',col=1)
#
#theta=0.8;
#res=dgeopoi(x,lambda=mean*theta, theta=theta);
#lines(0:x,res,t='l',col=2)
#
#theta=0.7;
#res=dgeopoi(x,lambda=mean*theta, theta=theta);
#lines(0:x,res,t='l',col=3)
#
#theta=0.6;
#res=dgeopoi(x,lambda=mean*theta, theta=theta);
#lines(0:x,res,t='l',col=4)
#
#theta=0.5;
#res=dgeopoi(x,lambda=mean*theta, theta=theta);
#lines(0:x,res,t='l',col=1)
#
#theta=0.4;
#res=dgeopoi(x,lambda=mean*theta, theta=theta);
#lines(0:x,res,t='l',col=2)
#
#theta=0.3;
#res=dgeopoi(x,lambda=mean*theta, theta=theta);
#lines(0:x,res,t='l',col=3)
#
#theta=0.2;
#res=dgeopoi(x,lambda=mean*theta, theta=theta);
#lines(0:x,res,t='l',col=4)
#
#theta=0.1;
#res=dgeopoi(x,lambda=mean*theta, theta=theta);
#lines(0:x,res,t='l',col=1)
