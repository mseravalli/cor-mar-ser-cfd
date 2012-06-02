#include "boundary_val.h"


/**
 * The boundary values of the problem are set.
 */
void boundaryvalues(
  int il,
  int ir,
  int jt,
  int jb,
  int omg_i,
  int omg_j,
  int iproc,
  int jproc,
  double **U,
  double **V
)
{
    int i = 0;
    int j = 0;

    /* set the boundaries only for the correct subdomains */

    /* left */
    if (omg_i == 0) {
        for (j = 1; j <= jt - jb; ++j) {
            U[0][j] = 0;
            V[0][j] = -V[1][j];
        }
    }

    /* right */
    if (omg_i == iproc - 1) {
        for (j = 1; j <= jt - jb; ++j) {
            U[ir][j] = 0;
            V[ir + 1][j] = -V[ir][j];
        }
    }

    /* bottom */
    if (omg_j == 0) {
        for (i = 1; i <= ir - il; ++i) {
            U[i][0] = -U[i][1];
            V[i][0] = 0;
        }
    }

    /* top */
    if (omg_j == jproc - 1) {
        for (i = 1; i <= ir - il; ++i) {
            U[i][jt + 1] = 2 - U[i][jt];
            V[i][jt] = 0;
        }
    }
}
