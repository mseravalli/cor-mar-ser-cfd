#include "boundary_val.h"


/**
 * The boundary values of the problem are set.
 */
void boundaryvalues(
  int imax,
  int jmax,
  double **U,
  double **V
)
{
    int i = 0;
    int j = 0;

    /** formula (14) **/
    for(j = 1; j <= jmax; ++j)
    {
        U[0][j]    = 0;
        U[imax+1][j] = 0;
    }

    for(i = 1; i <= imax; ++i)
    {
        U[i][0]    = 0;
        U[i][jmax+1] = 0;
    }

    /** formula (15) **/
    for(j = 1; j <= jmax; ++j)
    {
        V[0][j]      = -1.0 * V[1][j];
        V[imax+1][j] = -1.0 * V[imax][j];
    }

    for(i = 1; i <= imax; ++i)
    {
        V[i][0]      = -1.0 * V[i][1];
        V[i][jmax+1] = -1.0 * V[i][jmax];
    }

}
