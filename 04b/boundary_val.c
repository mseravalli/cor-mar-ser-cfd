#include "boundary_val.h"
#include "mpi.h"


/**
 * The boundary values of the problem are set.
 */
void boundaryvalues(
  int imax,
  int jmax,
  int rank_l,
  int rank_r,
  int rank_b,
  int rank_t,
  double **U,
  double **V
)
{
    int i = 0;
    int j = 0;

    /**** set the left boundary ****/
    if (rank_l == MPI_PROC_NULL) {
        for (j = 0; j <= jmax + 1; ++j) {
            U[1][j] = 0;
        }
        for (j = 0; j <= jmax + 2; ++j) {
            V[0][j] = -V[1][j];
        }
    }

    /**** set the right ****/
    if (rank_r == MPI_PROC_NULL) {
        for (j = 0; j <= jmax + 1; ++j) {
            U[imax+1][j] = 0;
        }
        for (j = 0; j <= jmax + 2; ++j) {
            V[imax+1][j] = -V[imax][j];
        }
    }

    /**** set the bottom boundary ****/
    if (rank_b == MPI_PROC_NULL) {
        for (i = 0; i <= imax + 2; ++i){
            U[i][0] = -U[i][1];
        }
        for (i = 0; i <= imax + 1; ++i) {
            V[i][1] = 0;
        }
    }

    /**** set the top boundary ****/
    if (rank_t == MPI_PROC_NULL) {
        for (i = 0; i <= imax + 2; ++i){
            U[i][jmax+1] = 2-U[i][jmax];
        }
        for (i = 0; i <= imax + 1; ++i) {
            V[i][jmax+1] = 0;
        }
    }

}
