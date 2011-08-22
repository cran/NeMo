
# Generate Minet graph miture
# alpha : proportion of nodes per class. The number of classes  is len(alpha)
# Pi : Pi[i,j] (symmetric/unsymmetric) matrix probabilities to connect nodes
#     from class i to nodes belonging to class j 
# symmetric : specifies if the returned edges list is undirected
# diag : if TRUE removes loops fron edge list
#
# Return an edge list[2, n_edges]
#
MixNet <- function(N, Alpha, Pi, storage="matrix", directed = FALSE, diag=FALSE){

  # Case where Pi is a scalar
  Pi <- as.matrix( Pi )
  

  testMixNetParameters( Alpha, Pi )
  
  # Test Pi matrix
  # Inv. ??? Juin 2010 Pi <- t(Pi)

  Alpha.n <- rowSums( rmultinom( N, 1, Alpha) )

  
  # n nodes is used to take into account the rounding effects
  n <- sum( Alpha.n)
    
  start.index <- vector()
  start.index[1] <- 1
  
  for( i in 1:( length(Alpha.n) )) {
    start.index[i+1] <- start.index[i] +  Alpha.n[i] 
  }
  adj <- matrix( 0, n, n)   
  for( i in 1:length(Alpha.n) ) {
    for( j in 1:length(Alpha.n) ) {
      i.beg = start.index[i] 
      i.end = start.index[i+1]-1
      j.beg = start.index[j] 
      j.end = start.index[j+1]-1
      if( (Alpha.n[i] != 0) & (Alpha.n[j] != 0) )
      adj[ i.beg:i.end, j.beg:j.end] <- matrix(
                             rbinom( Alpha.n[i] * Alpha.n[j], 1, Pi[i,j]),
                                                  Alpha.n[i],  Alpha.n[j]) 
    }
  }

  # 
  # Building Edge matrix
  #

  # Diagonnal terms
  if(!diag){
    for (i in 1:n){
    adj[i , i] = 0   
    }
  }
  if (storage == "edges") {
    if ( !directed  ){
      m <- t( which( (adj==1) & (upper.tri(adj)), arr.ind=TRUE) )
    } else {
      m <- t( which( (adj==1) , arr.ind=TRUE) )
    }
  } else if( storage == "matrix" ) {
    if ( !directed ) {
      m <- array( 0, dim( adj ) ) 
      adj[ lower.tri(adj) ] <- 0
      m <- adj + t( adj )
    } else {
      m <- adj
    }
  }
  return(m)
}

# Generate an EDD graph miture
# DegreeDist : Degree Distribution (matrix/vector) .
#              The indexes are shifted to used the 0-degree
#              DegreeDist[1] corresponds to zero degree value.
#              and so on.
# symmetric : specifies if the returned edges list is undirected
# diag : if TRUE removes loops fron edge list
#
# Return an edge list or an adjacency matrix
#
EDD <- function(N, deg.dist, storage="matrix", directed = FALSE, diag=FALSE){

  testEDDParameters( deg.dist, directed )
  
  if (directed) {
    n.degree <-  dim(deg.dist)[1]
  } else {
    n.degree <- length( deg.dist )
    # Duplicate distribution for Incoming and and
    # Outcoming arcs.
    deg.dist <- matrix( deg.dist, n.degree, 2 )
  }
    
  
  # Distribution Normalization
  sum.degree <-  colSums( deg.dist[,] )
  
  if (sum( sum.degree) == 0) {
    cat("Bad degree distribution deg.dist \n"  )
    return(NULL)
  }
  
  deg.dist[,1] <- deg.dist[,1]/sum.degree[1]
  deg.dist[,2] <- deg.dist[,2]/sum.degree[2]

  DegreeMean <- vector()
  DegreeMean[1] <- sum ( deg.dist[ ,1] %*% ( 0:(n.degree-1) ))
  DegreeMean[2] <- sum ( deg.dist[ ,2] %*% ( 0:(n.degree-1) ))
  gamma.array <-  1.0 / (DegreeMean * ( N - 1) );

  if( DegreeMean[1] != DegreeMean[2] ) {
    cat( "Bad distribution : expected degree in != expected degree out \n" )
    return(NULL)
  }

  # Search Degree max
  i <- dim( deg.dist ) [1]
  while ( deg.dist[i, 1] == 0 ) {
    i <- i - 1
  }
  max.deg.in <- i

  i <- dim( deg.dist ) [1]
  while ( deg.dist[i, 2] == 0 ) {
    i <- i - 1
  }
  max.deg.out <- i

  
  # The proba must be positive
  # Xij = B( k.DiDj ),
  #   with k = min{ 1/[(n-1) E(D)], 1/( max.deg.out * max.deg.in ) )
  gamma <- min ( gamma.array[1], 1/( max.deg.out * max.deg.in ) )

  #
  # Generate node degrees (In and Out) 
  #
  
  # In
  deg.trial <- rmultinom( N, 1, deg.dist[ ,1]) 
  # Take only the index which corresponds to the degre
  # and substract 1
  DegIn <- ( which ( deg.trial == 1 , arr.ind=TRUE ) [,1] ) - 1

  # Out
  if ( directed ) {
    deg.trial <- rmultinom( N, 1, deg.dist[ , 2 ]) 
    # Take only the index which corresponds to the degre
    # and substract 1
    DegOut <- which ( deg.trial == 1 , arr.ind=TRUE ) [,1]  - 1

  } else {
    # To have the same distribution for the undirected case 
    DegOut <- DegIn
  }


  # Compute the Adjacency matrix 
  adj <- matrix( 0, N, N)

  for( i in 1:N ) {
    VA <- runif( N )
    proba <- gamma * DegIn[i] * DegOut[]
     adj[ i, ] = as.integer( VA < proba )
  }
  
  # 
  # Building an Adjacency matrix or edge list
  #

  # Diagonnal terms
  if(!diag){
    for (i in 1:N){
    adj[i , i] <- 0   
    }
  }
  if (storage == "edges") {
    if ( !directed ){
      m <- t( which( (adj==1) & (upper.tri(adj)), arr.ind=TRUE) )
    } else {
      m <- t( which( (adj==1) , arr.ind=TRUE) )
    }
  } else if (storage == "matrix") {
    if ( !directed ) {
      m <- array( 0, dim( adj ) ) 
      adj[ lower.tri(adj) ] <- 0
      m <- adj + t( adj )
    } else {
      m <- adj
    }
  }
  return(m)
}

AdjMat2Edges<-function(x, directed=FALSE, verbose=FALSE ) {

  # Warning : return only connected nodes
  if (dim(x)[1] == dim(x)[2]){
    
    # Adjacency matrix case
    nbrNodes<-dim(x)[1]
    NotConnected <- which( (rowSums( x ) + colSums(x)) == 0)
    if ( length(NotConnected) > 0 ) {
      # GG inv
      # Not necessary x <- x[ - NotConnected, - NotConnected ]
       if (verbose)
         cat("Some nodes are not connected to the network \n")
    } 
    if ( directed ) {
      m <- t( which( (t(x)==1) , arr.ind=TRUE) )
    } else {
      m <- t( which( (x==1) & (upper.tri(x)), arr.ind=TRUE) )
    }
  } else {
    cat( "Error: Bad matrix dimensions \n" )
    stop()
  }

  return( m )
}

# ????
display <- function(motif, directed=FALSE ) {

  par(mar=c(0, 4, 1, 1))
  par(mfrow=c(4, 2))

   if ( ! is.list( motif ) )
    motif <- list( motif )
   for( m in 1:length(motif) ) {
     if( directed ) {
       gplot( motif[[m]], gmode='digraph', label=0:( dim(motif[[m]])[1]-1 ),
             edge.col = 'blue', arrowhead.cex = 3,
             main= paste( "motif", m) )
     } else {
       gplot( motif[[m]],  gmode='graph', label=0:( dim(motif[[m]])[1]-1 ),
             edge.col = 'blue',
             main= paste( "motif", m) )
     }       
   }
}

