#include "boundary_val.h"
#include "constants.h"
#include <string.h>

/**
 * The boundary values of the problem are set.
 */
void boundaryvalues(
  int imax,
  int jmax,
  int kmax,
  double dx,
  double dy,
  int wl,
  int wr,
  int wt,
  int wb,
  double** U,
  double** V,
  double** F,
  double** G,
  double** P,
  double*** C,
  int** Flag
)
{
    int i = 0;
    int j = 0;
    int k = 0;

    /* probably those values should be passed as parameters */
    /* MUST be Dynamic */
    /* stored as left, right, bottom, top */
    double C_bound[4][4] = {{0.0, 0.0, 0.0, 1.0},
                            {0.0, 0.0, 1.0, 0.0},
                            {0.0, 1.0, 0.0, 0.0},
                            {1.0, 0.0, 0.0, 0.0}};

    /******** EXTERNAL WALLS START ********/

    /** left wall **/
    if(wl==1) /** no slip **/
    {
         for(j = 1; j <= jmax; ++j)
         {
             U[0][j] = 0;
             V[0][j] = -1.0*V[1][j];
             for (k = 0; k < kmax; k++) {
                 C[k][0][j] = 2*C_bound[k][0]*((j-0.5)*dy) - C[k][1][j];
             }
         }
    }
    else if (wl==2) /** free slip **/
    {
         for(j = 1; j <= jmax; ++j)
         {
             U[0][j] = 0;
             V[0][j] = V[1][j];
             for (k = 0; k < kmax; ++k) {
                 C[k][0][j] = C[k][1][j] + dx*C_bound[k][0]*((j-0.5)*dy);
             }
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
             for (k = 0; k < kmax; k++) {
                 C[k][imax+1][j] = 2*C_bound[k][1]*((j-0.5)*dy) - C[k][imax][j];
             }
         }
    }
    else if (wr==2) /** free slip **/
    {
         for(j = 1; j <= jmax; ++j)
         {
             U[imax][j] = 0;
             V[imax+1][j] = V[imax][j];
             for (k = 0; k < kmax; ++k) {
                 C[k][imax+1][j] = C[k][imax][j] + dx*C_bound[k][1]*((j-0.5)*dy);
             }
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

    /** lower wall **/
    if(wb==1) /** no slip **/
    {
         for(i = 1; i <= imax; ++i)
         {
             V[i][0] = 0;
             U[i][0] = -1.0*U[i][1];
             for (k = 0; k < kmax; k++) {
                 C[k][i][0] = 2*C_bound[k][2]*((i-0.5)*dx) - C[k][i][1];
             }
         }
    }
    else if (wb==2) /** free slip **/
    {
         for(i = 1; i <= imax; ++i)
         {
             V[i][0] = 0;
             U[i][0] = U[i][1];
             for (k = 0; k < kmax; ++k) {
                 C[k][i][0] = C[k][i][1] + dy*C_bound[k][2]*((i-0.5)*dx);
             }
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

    /** upper wall **/
    if(wt==1) /** no slip **/
    {
         for(i = 1; i <= imax; ++i)
         {
             V[i][jmax] = 0;
             U[i][jmax+1] = -1.0*U[i][jmax];
             for (k = 0; k < kmax; k++) {
                 C[k][i][jmax+1] = 2*C_bound[k][3]*((i-0.5)*dx) - C[k][i][jmax];
             }
         }
    }
    else if (wt==2) /** free slip **/
    {
         for(i = 1; i <= imax; ++i)
         {
             V[i][jmax] = 0;
             U[i][jmax+1] = U[i][jmax];
             for (k = 0; k < kmax; ++k) {
                 C[k][i][jmax+1] = C[k][i][jmax] + dy*C_bound[k][3]*((i-0.5)*dx);
             }
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
    
    /******** EXTERNAL WALLS END ********/


    /******** INSTERNAL OBJECTS START ********/

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
                F[i-1][j]=U[i-1][j];
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

    /******** INSTERNAL OBJECTS END ********/

 
}

void spec_boundary_val(
    char *problem,
    int imax,
    int jmax,
    double dx,
    double dy,
    double Re,
    double deltaP,
    double **U,
    double **V,
    double **P
    )
{
    int j = 0;
    /*double dpdx;*/

    if (strcmp(problem, "karman") == 0) {
        for(j = 1; j <= jmax; ++j)
        {
            U[0][j] = 1.0;
            V[0][j] = -V[1][j];
        }
    }

}

