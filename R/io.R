readAdjacency <- function( file, directed=FALSE, diag=TRUE, verbose=FALSE ) {
# - directed=TRUE means that the read matrix is not transformed.
#   An arc is built from the first field to the second field 
# - directed=FALSE means that the matrix is symmetrized
#   duplicated edges are removed
# - edge indexes are stored in decreasing  order
# - "space" is the field separator

  nodenames <- vector()
  n.nodes   <- 0
  edges     <- array() 
  n.edges   <- 0
  src       <- vector()
  dst       <- vector()
  
  if (verbose) cat("\nReading file ...",  "\n")
  
  lines <- readLines(file)
  for ( i in 1:length( lines ) ) {
    fields <- unlist(strsplit(lines[i],"[[:space:]]+",perl=TRUE))

    # Add in the node list if a new one
    if ( !( fields[1] %in% nodenames ) ) {
      n.nodes <- n.nodes+1
      nodenames[n.nodes] <-  fields[1]  
    }

    if ( !( fields[2] %in% nodenames ) ) {
      n.nodes <- n.nodes+1
      nodenames[n.nodes] <-  fields[2]  
    }

    # Add an edge
    src.test <- which( nodenames == fields[1] )[1] 
    dst.test <- which( nodenames == fields[2] )[1]

    # Test edge is already in list
    # ----------------------------

    # Diagonal terms
    if( diag | (src.test != dst.test) ) {
      if( length( which( dst[ which( src == src.test) ] == dst.test ) ) == 0) {
        
        if ( ! directed ) {
          if( length( which( dst[ which( src == dst.test) ] == src.test ) ) == 0) {
            n.edges <- n.edges + 1
            if ( src.test >= dst.test ) {
              src[ n.edges] <- src.test
              dst[ n.edges] <- dst.test     
            } else {
              src[ n.edges] <- dst.test
              dst[ n.edges] <- src.test
            }
          } else {
            if( verbose )
              cat( "Warning: remove line", i," ->", lines[i], "\n") 
          }
        } else {
          n.edges <- n.edges + 1
          src[ n.edges] <- src.test
          dst[ n.edges] <- dst.test
        }
      } else {
        if( verbose )
          cat( "Warning: remove line", i," ->", lines[i], "\n") 
      }
    } else {
      if( verbose )
        cat( "Warning: remove line", i," ->", lines[i], "\n") 

    }
  }
  invisible( list( directed=directed, n.nodes = n.nodes, nodenames = nodenames,
               n.edges =n.edges, edges = rbind(src, dst) ) )   
}
 
getDegrees <- function( x, in.out=FALSE ) {
  ## ??? S3 method

  Degrees.in <- rep( 0, x$n.nodes )
  Degrees.out <- rep( 0, x$n.nodes )
  for ( i in 1:x$n.edges ) {
    node = x$edges[1, i]
    Degrees.out[ node ] = Degrees.out[ node ] + 1
    node = x$edges[2, i]
    Degrees.in[ node ] = Degrees.in[ node ] + 1
  }
  res <- Degrees.out +  Degrees.in 

  if( in.out ) {
    res <- list( d.in= Degrees.in, d.out= Degrees.out, d.sum=res )
  } 
  return( res )
  
}

getDist <- function( x ) {
# dist[1] <- 0 degree
# dist[2] <- 1 degree
# ...
# dist[max(x)+1] <- max(x)
  
  x <- as.vector( x )
  n <- length( x )
  x.max <- max( x )
  dist <- rep( 0, x.max + 1)
  for ( i in 1:n ) {
    dist[ x[i]+1 ] <-  dist[ x[i]+1 ] + 1
  }     
  invisible( dist )
}
         

Edges2AdjMat<-function(x, symmetrize=FALSE ) {
  # ??? x object
  n <- max( x )

  mat <- matrix(0, n, n)
  mat[ cbind( x[1,], x[2,] ) ] <- 1
  if (symmetrize)
    mat[ cbind( x[2,], x[1,] ) ] <- 1
  invisible( mat )
}
