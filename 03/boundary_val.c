#include "boundary_val.h"
#include "constants.h"
#include <string.h>

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
  double **V,
  double **F,
  double **G,
  int** Flag
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


    /* set the values of the inner obstacles */
    for (i = 1; i <= imax; ++i) { 
        for (j = 1; j <= jmax; ++j) {

            if (Flag[i][j] == B_N)
            {
                V[i][j]=0;
                U[i-1][j]=-1.0*U[i-1][j+1];
                U[i][j]=-1.0*U[i][j+1];
                G[i][j]=V[i][j];
            }
            
            else if (Flag[i][j] == B_W)
            {
                U[i-1][j]=0;
                V[i][j]=-1.0*V[i-1][j];
                V[i][j-1]=-1.0*V[i-1][j-1];
                F[i-1][j]=U[i-1][j];
            }

            else if (Flag[i][j] == B_S)
            {
                V[i][j-1]=0;
                U[i-1][j]=-1.0*U[i-1][j-1];
                U[i][j]=-1.0*U[i][j-1];
                G[i][j-1]=V[i][j-1];
            }

            else if (Flag[i][j] == B_O)
            {
                U[i][j]=0;
                V[i][j]=-1.0*V[i+1][j];
                V[i][j-1]=-1.0*V[i+1][j-1];
                F[i][j]=U[i][j];
            }

            else if (Flag[i][j] == B_NO) 
            {
                U[i][j]=0;
                V[i][j]=0;
                U[i-1][j]=-1.0*U[i-1][j+1];
                V[i][j-1]=-1.0*V[i+1][j-1];
                F[i][j]=U[i][j];
                G[i][j]=V[i][j];
            }

            else if (Flag[i][j] == B_NW) 
            {
                U[i-1][j]=0;
                V[i][j]=0;
                U[i][j]=-1.0*U[i][j+1];
                V[i][j-1]=-1.0*V[i-1][j-1];
                F[i-1][j]=U[i][j];
                G[i][j]=V[i][j];
            }

            else if (Flag[i][j] == B_SO) 
            {
                U[i][j]=0;
                V[i][j-1]=0;
                U[i-1][j]=-1.0*U[i-1][j-1];
                V[i][j]=-1.0*V[i+1][j];
                F[i][j]=U[i][j];
                G[i][j-1]=V[i][j-1];
            }

            else if (Flag[i][j] == B_SW) 
            {
                U[i-1][j]=0;
                V[i][j-1]=0;
                U[i][j]=-1.0*U[i][j-1];
                V[i][j]=-1.0*V[i-1][j];
                F[i-1][j]=U[i-1][j];
                G[i][j-1]=V[i][j-1];
            }
        }
    }

 
}

void spec_boundary_val(
    char *problem,
    int imax,
    int jmax,
    double dx,
    double dy,
    double Re,
    double **U,
    double **V,
    double **P
    )
{
    int j = 0;
    double dpdx;

    if (strcmp(problem, "karman") == 0) {
        for(j = 1; j <= jmax; ++j)
        {
            U[0][j] = 1.0;
            V[0][j] = 0.0;
        }
    } else if (strcmp(problem, "plane") == 0) {
        for(j = 1; j <= jmax; ++j)
        {
            dpdx = (P[1][j] - P[0][j]) / dx;
            U[0][j] = -0.5 * Re * (j - 1)*dy * (j - jmax)*dy * dpdx;
            V[0][j] = 0.0;
        }
    }

}

