#include "boundary_val.h"


/**
 * The boundary values of the problem are set.
 */
void boundaryvalues(
  int imax,
  int jmax,
  int wl,
  int wr,
  int wt,
  int wb,
  double **U,
  double **V
)
{
    int i = 0;
    int j = 0;

    /** left wall **/
    if(wl==1) /** no slip **/
    {
         for(j = 1; j <= jmax; ++j)
         {
             U[0][j] = 0;
             V[0][j] = -1.0*V[1][j];
         }
    }
    else if (wl==2) /** free slip **/
    {
         for(j = 1; j <= jmax; ++j)
         {
             U[0][j] = 0;
             V[0][j] = V[1][j];
         }
    }
    else if (wl==3) /** Outflow **/
    {
         for(j = 1; j <= jmax; ++j)
         {
             U[0][j] = U[1][j];
             V[0][j] = V[1][j];
         }
    }

    /** right wall **/
    if(wr==1) /** no slip **/
    {
         for(j = 1; j <= jmax; ++j)
         {
             U[imax][j] = 0;
             V[imax+1][j] = -1.0*V[imax][j];
         }
    }
    else if (wr==2) /** free slip **/
    {
         for(j = 1; j <= jmax; ++j)
         {
             U[imax][j] = 0;
             V[imax+1][j] = V[imax][j];
         }
    }
    else if (wr==3) /** Outflow **/
    {
         for(j = 1; j <= jmax; ++j)
         {
             U[imax][j] = U[imax-1][j];
             V[imax+1][j] = V[imax][j];
         }
    }

    /** upper wall **/
    if(wt==1) /** no slip **/
    {
         for(i = 1; i <= imax; ++i)
         {
             V[i][jmax] = 0;
             U[i][jmax+1] = -1.0*U[i][jmax];
         }
    }
    else if (wt==2) /** free slip **/
    {
         for(i = 1; i <= imax; ++i)
         {
             V[i][jmax] = 0;
             U[i][jmax+1] = U[i][jmax];
         }
    }
    else if (wt==3) /** Outflow **/
    {
         for(i = 1; i <= imax; ++i)
         {
             U[i][jmax+1] = U[i][jmax];
             V[i][jmax] = V[i][jmax-1];
         }
    }


    /** lower wall **/
    if(wb==1) /** no slip **/
    {
         for(i = 1; i <= imax; ++i)
         {
             V[i][0] = 0;
             U[i][0] = -1.0*U[i][1];
         }
    }
    else if (wb==2) /** free slip **/
    {
         for(i = 1; i <= imax; ++i)
         {
             V[i][0] = 0;
             U[i][0] = U[i][1];
         }
    }
    else if (wb==3) /** Outflow **/
    {
         for(i = 1; i <= imax; ++i)
         {
             U[i][0] = U[i][1];
             V[i][0] = V[i][1];
         }
    }
 
}

void spec_boundary_val(
    char *problem,
    int imax,
    int jmax,
    double **U,
    double **V
    )
{
 /** TODO **/

}
