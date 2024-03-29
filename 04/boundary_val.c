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
        for (j = 0; j <= jt - jb + 2; ++j) {
            U[1][j] = 0;
            V[0][j] = -V[1][j];
        }
    }

    /* right */
    if (omg_i == iproc - 1) {
        for (j = 0; j <= jt - jb + 2; ++j) {
            U[ir-il+2][j] = 0;
            V[ir-il+2][j] = -V[ir-il+1][j];
        }
    }

    /* bottom */
    if (omg_j == 0) {
        for (i = 0; i <= ir - il + 2; ++i) {
            U[i][0] = -U[i][1];
            V[i][1] = 0;
        }
    }

    /* top */
    if (omg_j == jproc - 1) {
        for (i = 0; i <= ir - il + 2; ++i) {
            U[i][jt-jb+2] = 2 - U[i][jt-jb+1];
            V[i][jt-jb+2] = 0;
        }
    }
}
