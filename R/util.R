sub.param <- function (param, default=NULL, ... ) {

  if (missing(param)) { return(default) }

  l <- list(...)
  res <- l[[  which( names(l) %in% param)[1] ]] 

  if (is.null(res)) { res <- default }

  return ( res )

}


testMixNetParameters <- function ( Alpha, Pi ) {
  
  if(  is.null( Alpha ) | is.null( Pi ) ) {
    cat( "Error: for 'MixNet' model, Alpha and Pi are required \n")
    stop()
  } else if ( ! is.vector( Alpha ) | ! is.matrix( Pi ) ){
    cat( "Error: Alpha must be a vector and Pi must be a matrix \n")
    stop()
  } else if( ! all ( dim( Alpha ) == dim( Pi ) ) ) {
    cat( "Error: Alpha and Pi have not the same dimensions \n")
    stop()
  }
}

testEDDParameters <- function ( deg.dist, directed ) {
  
  if ( ! is.null( deg.dist ) ) {
    if( directed ) {
      if( is.matrix( deg.dist ) ) {
        if ( dim( deg.dist )[2] != 2 ) {
          cat( "Error: 'deg.dist' must be a 2 column matrix \n")
          stop()
        }
      } else {
        cat( "Error: 'deg.dist' must be a 2 column matrix \n")
        stop()
      }
    }
  }
}

testERParameters <- function ( Pi ) {
  
  if ( ! is.null( Pi ) ) {
    if ( ! is.vector( Pi ) ) {
      cat( "Error: Pi must be scalar for ER model\n")
      stop()
    } else if ( length( Pi ) > 1 ) {
      cat( "Error: Pi must be scalar for ER model\n")
      stop()
    }
  }  
}


display <- function(motif, directed=FALSE ) {

  par(mar=c(0, 4, 1, 1))
  par(mfrow=c(4, 2))

  if ( ! is.list( motif ) )
    motif <- list( motif )
  for( m in 1:length(motif) ) {
    if( directed ) {
       # label=0:( dim(motif[[m]])[1]-1 ),
      Gplot( motif[[m]], directed=directed, main= paste( "motif", m) )
    } else {
      Gplot( motif[[m]], directed=directed, main= paste( "motif", m) )
    }       
  }
}

# Set the path of the motif data base (C++ part)
setMotifsPath <- function() {
  motif.path = system.file("networks",package="NeMo")
  # cat( "motif.path = ", motif.path, "\n")
         .Call( "setMotifPath_C", motif.path,
                PACKAGE="NeMo")

}
