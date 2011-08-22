#include <iostream>
#include <string.h>
#include <stdint.h>
// Must be before R includes - conflict with length function
# include "graphical-models.h"
# include "motif-search.h"
# include <R.h>
# include <Rdefines.h>



using namespace std;

extern "C" {

SEXP  setMotifPath_C(  SEXP S_path ) {

  PROTECT(S_path = AS_CHARACTER(S_path));
  char *path = new char[ strlen(CHAR(STRING_ELT(S_path, 0))) ];
  strcpy(path, CHAR(STRING_ELT(S_path, 0)));
  
  setMotifDirectory( path );
  // cout << "C Path: " << path << endl;

  UNPROTECT(1);
  return(R_NilValue);
}

SEXP countNonExactMotif_C( SEXP S_G, SEXP S_G_nnodes, 
			   SEXP S_M, SEXP S_M_nnodes,
			   SEXP S_Directed ) {
  // S_G edge list of the network G
  // S_G_nnodes number of nodes of G
  // S_M edge list of the motif 
  // S_M_nnodes number of nodes of M
  // S_Directed non-zero integer if the graph is directed

  int n_prot = 0;

  // Arguments
  PROTECT( S_G        = AS_INTEGER(S_G)    ); n_prot++;
  PROTECT( S_G_nnodes = AS_INTEGER(S_G_nnodes) ); n_prot++;
  PROTECT( S_M        = AS_INTEGER(S_M)    ); n_prot++;
  PROTECT( S_M_nnodes = AS_INTEGER(S_M_nnodes) ); n_prot++;
  PROTECT( S_Directed = AS_INTEGER(S_Directed) ); n_prot++;

  // Passing argument by addresses
  int *G_edges = INTEGER_POINTER( S_G );
  int G_nnodes = INTEGER_POINTER( S_G_nnodes )[0];
  int *M_edges = INTEGER_POINTER( S_M );
  int M_nnodes = INTEGER_POINTER( S_M_nnodes )[0];
  bool Directed = ( INTEGER_POINTER( S_Directed )[0] != 0 );

  SparseAdjacency G( G_edges, G_nnodes, ! Directed );
  SparseAdjacency M( M_edges, M_nnodes, ! Directed );

  // FindMotif find( CountMode, &cout );
  FindMotif find( CountMode, 0 );

  find.find( G, M );

  int count = find.getCount();

  // Return List
  SEXP S_ReturnVal;


  // Allocate and store result
  PROTECT( S_ReturnVal  = allocVector( INTSXP, 1 )); n_prot++;
  INTEGER_POINTER( S_ReturnVal ) [0] = count;


  UNPROTECT(n_prot);

  return(S_ReturnVal);

}

SEXP countAllNonExactMotif_C( SEXP S_G, SEXP S_G_nnodes, 
			      SEXP S_M_nnodes,
			      SEXP S_Directed ) {
  // S_G edge list of the network G
  // S_G_nnodes number of nodes of G
  // S_M edge list of the motif 
  // S_M_nnodes number of nodes of M
  // S_Directed non-zero integer if the graph is directed

  int n_prot = 0;

  // Arguments
  PROTECT( S_G        = AS_INTEGER(S_G)    ); n_prot++;
  PROTECT( S_G_nnodes = AS_INTEGER(S_G_nnodes) ); n_prot++;
  PROTECT( S_M_nnodes = AS_INTEGER(S_M_nnodes) ); n_prot++;
  PROTECT( S_Directed = AS_INTEGER(S_Directed) ); n_prot++;

  // Passing argument by addresses
  int *G_edges = INTEGER_POINTER( S_G );
  int G_nnodes = INTEGER_POINTER( S_G_nnodes )[0];
  int M_nnodes = INTEGER_POINTER( S_M_nnodes )[0];
  bool Directed = ( INTEGER_POINTER( S_Directed )[0] != 0 );

  SparseAdjacency G( G_edges, G_nnodes, ! Directed );

  FindMotif find(CountMode, 0);
  vector <int*> motifs;
  vector <int64_t > counts;
  int k = M_nnodes;
  buildAndfindAllMotif( G, k, find, motifs, counts, ReadMotifs, NonExactMode );

  SEXP S_ReturnList;
  PROTECT( S_ReturnList=allocVector(VECSXP, 2 * counts.size() ) ); n_prot++;
  int i, j;


  // Index list
  int nbr_motif =  motifs.size();
  int m = 0;
  for( vector<int*>::iterator it_m = motifs.begin() ;  
       it_m != motifs.end(); m++, it_m++ ) {
    int* index = *it_m;
    SEXP S_Adj;
    SEXP S_Counts;
    PROTECT( S_Counts = allocVector( REALSXP, 1 ) ); n_prot++; 
    double *Counts = NUMERIC_POINTER( S_Counts );

    PROTECT( S_Adj = allocVector( INTSXP, k*k ) ); n_prot++;
    int *adj = INTEGER_POINTER( S_Adj );

    // Set to 0 Adj matrix
    for( int *it = adj, *it_end = &adj[k*k]; it != it_end; *it++ = 0 ) ;

    for( int *it = index ; *it != NOT_INDEX; ) {
      i = *it++;
      j = *it++;
      adj[ i + j * k] = 1;
      if( ! Directed ) {
	adj[ j + i * k] = 1;
      }
    }
    delete [] index;

    Counts[0] = counts[m]; 
    SET_VECTOR_ELT( S_ReturnList, m, S_Adj );
    SET_VECTOR_ELT( S_ReturnList, m + nbr_motif, S_Counts );
  }

  UNPROTECT( n_prot );

  return( S_ReturnList );

}


SEXP computeERMoment_C( SEXP S_N,
			SEXP S_Pi, 
			SEXP S_M, SEXP S_M_nedges,
			SEXP S_Directed ) {
  // S_N Number of nodes. 
  // S_Pi connectivity probability.
  // S_M adjacency matrix of the motif 
  // S_M_nedges edge number in the motif
  // S_Directed non-zero integer if the graph is directed

  int n_prot = 0;

  // Arguments
  PROTECT( S_N        = AS_INTEGER(S_N)        ); n_prot++;
  PROTECT( S_Pi       = AS_NUMERIC(S_Pi)       ); n_prot++;
  PROTECT( S_M        = AS_INTEGER(S_M)        ); n_prot++;
  PROTECT( S_M_nedges = AS_INTEGER(S_M_nedges) ); n_prot++;
  PROTECT( S_Directed = AS_INTEGER(S_Directed) ); n_prot++;

  int64_t N     = INTEGER_POINTER( S_N )[0];
  double Pi       = NUMERIC_POINTER( S_Pi )[0];
  int *M          = INTEGER_POINTER( S_M );
  size_t M_nedges = INTEGER_POINTER( S_M_nedges )[0];
  bool Directed   = ( INTEGER_POINTER( S_Directed )[0] != 0 );

  // Allocate storage for results
  SEXP S_Res;
  PROTECT( S_Res = NEW_NUMERIC( 2 ) ); n_prot++;
  double *Res = NUMERIC_POINTER( S_Res );

  // Transform to int
  int max = 0; 
  int *M_t = new int[M_nedges*2+2];
  for ( size_t i=0; i < M_nedges*2; i++) {
    M_t [i] = M[i];
    max = MAX( max,  M_t[i] );
  }
  M_t [M_nedges*2]   = NOT_INDEX;
  M_t [M_nedges*2+1] = NOT_INDEX;
  int nnodes = max + 1;

  SparseAdjacency motif( M_t, nnodes, !Directed );

  NonRedondantPermutation NRPMotif( motif );

  ErdosRenyi ER(N, Pi, Directed);

  Res[0] = ER.computeMean( motif, NRPMotif );
  Res[1] = ER.computeMoment2( motif, NRPMotif );

  UNPROTECT(n_prot);

  delete [] M_t;

  return(S_Res);
}

SEXP computeEDDMoment_C( SEXP S_N, 
			    SEXP S_DegreeDist, SEXP S_NDegree,
			    SEXP S_M, SEXP S_M_nedges,
			    SEXP S_Directed ) {
  // S_N Number of nodes. 
  // S_DegreeDist
  // S_NDegrees
  // S_M adjacency matrix of the motif 
  // S_M_nedges edge number in the motif
  // S_Directed non-zero integer if the graph is directed

  int n_prot = 0;

  // Arguments
  PROTECT( S_N          = AS_INTEGER(S_N)         ); n_prot++;
  PROTECT( S_DegreeDist = AS_NUMERIC(S_DegreeDist)); n_prot++;
  PROTECT( S_NDegree    = AS_INTEGER(S_NDegree)   ); n_prot++;
  PROTECT( S_M          = AS_INTEGER(S_M)         ); n_prot++;
  PROTECT( S_M_nedges   = AS_INTEGER(S_M_nedges)  ); n_prot++;
  PROTECT( S_Directed   = AS_INTEGER(S_Directed)  ); n_prot++;

  int64_t N        = INTEGER_POINTER( S_N )[0];
  int    NDegree     = INTEGER_POINTER( S_NDegree )[0];
  double *DegreeDist = NUMERIC_POINTER( S_DegreeDist );
  int *M             = INTEGER_POINTER( S_M );
  size_t M_nedges    = INTEGER_POINTER( S_M_nedges )[0];
  bool Directed      = ( INTEGER_POINTER( S_Directed )[0] != 0 );

  // Allocate storage for results
  SEXP S_Res;
  PROTECT( S_Res = NEW_NUMERIC( 2 ) ); n_prot++;
  double *Res = NUMERIC_POINTER( S_Res );

  // Transform to size_t
  int max = 0; 
  int *M_t = new int[M_nedges*2+2];
  for ( size_t i=0; i < M_nedges*2; i++) {
    M_t [i] = M[i];
    max = MAX( max,  M_t[i] );
  }
  M_t [M_nedges*2]   = NOT_INDEX;
  M_t [M_nedges*2+1] = NOT_INDEX;
  int nnodes = max + 1;

  SparseAdjacency motif( M_t, nnodes, !Directed );

  NonRedondantPermutation NRPMotif( motif );

  EDD edd(N, DegreeDist, NDegree, Directed);

  Res[0] = edd.computeMean( motif, NRPMotif );
  Res[1] = edd.computeMoment2( motif, NRPMotif );

  UNPROTECT(n_prot);

  delete [] M_t;

  return(S_Res);
}

SEXP computeMixNetMoment_C( SEXP S_N, SEXP S_nclass,
			    SEXP S_Alpha, SEXP S_Pi, 
			    SEXP S_M, SEXP S_M_nedges,
			    SEXP S_Directed ) {
  // S_N Number of nodes. 
  // S_nclass Number of classes. 
  // S_Alpha number of nodes per classes
  // S_Pi  Connectivity probabilities between classes.
  // S_M adjacency matrix of the motif 
  // S_M_nedges edge number in the motif
  // S_Directed non-zero integer if the graph is directed

  int n_prot = 0;

  // Arguments
  PROTECT( S_N        = AS_INTEGER(S_N)        ); n_prot++;
  PROTECT( S_nclass   = AS_INTEGER(S_nclass)   ); n_prot++;
  PROTECT( S_Alpha    = AS_NUMERIC(S_Alpha)    ); n_prot++;
  PROTECT( S_Pi       = AS_NUMERIC(S_Pi)       ); n_prot++;
  PROTECT( S_M        = AS_INTEGER(S_M)        ); n_prot++;
  PROTECT( S_M_nedges = AS_INTEGER(S_M_nedges) ); n_prot++;
  PROTECT( S_Directed = AS_INTEGER(S_Directed) ); n_prot++;

  int64_t N     = INTEGER_POINTER( S_N )[0];
  int    nclass   = INTEGER_POINTER( S_nclass )[0];
  double *Alpha   = NUMERIC_POINTER( S_Alpha );
  double *Pi      = NUMERIC_POINTER( S_Pi );
  int *M          = INTEGER_POINTER( S_M );
  size_t M_nedges = INTEGER_POINTER( S_M_nedges )[0];
  bool Directed   = ( INTEGER_POINTER( S_Directed )[0] != 0 );

  // Allocate storage for results
  SEXP S_Res;
  PROTECT( S_Res = NEW_NUMERIC( 2 ) ); n_prot++;
  double *Res = NUMERIC_POINTER( S_Res );

  // Transform to size_t
  int max = 0; 
  int *M_t = new int[M_nedges*2+2];
  for ( size_t i=0; i < M_nedges*2; i++) {
    M_t [i] = M[i];
   max = MAX( max,  M_t[i] );
  }
  M_t [M_nedges*2]   = NOT_INDEX;
  M_t [M_nedges*2+1] = NOT_INDEX;
  int nnodes = max + 1;

  SparseAdjacency motif( M_t, nnodes, !Directed );

  NonRedondantPermutation NRPMotif( motif );

  Mixnet mixnet(N, nclass, Alpha, Pi, Directed);

  Res[0] = mixnet.computeMean( motif, NRPMotif );
  Res[1] = mixnet.computeMoment2( motif, NRPMotif );

  UNPROTECT(n_prot);

  delete [] M_t;
  return(S_Res);

}

}
