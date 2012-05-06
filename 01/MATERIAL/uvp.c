#include "boundary_val.h"
#include "disc.h"
#include "helper.h"
#include "matrix_op.h"
#include "uvp.h"
#include <math.h>


void calculate_fg(
  double Re,
  double GX,
  double GY,
  double alpha,
  double dt,
  double dx,
  double dy,
  int imax,
  int jmax,
  double **U,
  double **V,
  double **F,
  double **G
)
{
    
    /******** VARIABLE DECLARATION START ********/

    double** _d2udx2;
    double** _d2udy2;
    double** _du2dx ;
    double** _duvdy ;

    double** _d2vdx2;
    double** _d2vdy2;
    double** _duvdx ;
    double** _dv2dy ;
    
    int i;
    int j;

    /******** VARIABLE DECLARATION END ********/


    /******** CALCULATE F START ********/
    
    boundaryvalues(imax, jmax, U, V);

    /******** Calculate the single derivatives ********/

    /* calculate derivatives */
    _d2udx2 = d2udx2(U, dx,   0, imax + 1, 0, jmax + 1);
    _d2udy2 = d2udy2(U, dy,   0, imax + 1, 0, jmax + 1);
    _du2dx  = du2dx(U, dx,    0, imax + 1, 0, jmax + 1);
    _duvdy  = duvdy(U, V, dy, 0, imax + 1, 0, jmax + 1);

    /******** Calculate F ********/
    add_mat(   _d2udx2, _d2udy2, 0, imax+1, 0, jmax+1, F);
    mult_scalar(F,       1/Re,   0, imax+1, 0, jmax+1, F);
    sub_mat(    F,       _du2dx, 0, imax+1, 0, jmax+1, F);
    sub_mat(    F,       _duvdy, 0, imax+1, 0, jmax+1, F);
    add_scalar( F,       GX,     0, imax+1, 0, jmax+1, F);
    mult_scalar(F,       dt,     0, imax+1, 0, jmax+1, F);
    add_mat(    F,       U,      0, imax+1, 0, jmax+1, F);

    /******** Free the derivatives ********/
    free_matrix(_d2udx2, 0, imax+1, 0, jmax+1);
    free_matrix(_d2udy2, 0, imax+1, 0, jmax+1);
    free_matrix(_du2dx,  0, imax+1, 0, jmax+1);
    free_matrix(_duvdy,  0, imax+1, 0, jmax+1);

    /******** CALCULATE F END ********/
    

    /******** CALCULATE G START ********/

    boundaryvalues(imax, jmax, U, V);

    /******** Calculate the single derivatives ********/

    /* calculate derivatives */
    _d2vdx2 = d2vdx2(U, dx,   0, imax+1, 0, jmax+1);
    _d2vdy2 = d2vdy2(U, dy,   0, imax+1, 0, jmax+1);
    _duvdx  = duvdx(U, V, dy, 0, imax+1, 0, jmax+1);
    _dv2dy  = dv2dy(U, dx,    0, imax+1, 0, jmax+1);

    /******** Calculate G ********/
    add_mat(   _d2vdx2, _d2vdy2, 0, imax+1, 0, jmax+1, G);
    mult_scalar(G,      1/Re,    0, imax+1, 0, jmax+1, G);
    sub_mat(    G,      _duvdx,  0, imax+1, 0, jmax+1, G);
    sub_mat(    G,      _dv2dy,  0, imax+1, 0, jmax+1, G);
    add_scalar( G,      GY,      0, imax+1, 0, jmax+1, G);
    mult_scalar(G,      dt,      0, imax+1, 0, jmax+1, G);
    add_mat(    G,      V,       0, imax+1, 0, jmax+1, G);

    /******** Free the derivatives ********/
    free_matrix(_d2vdx2,    0, imax+1, 0, jmax+1);
    free_matrix(_d2vdy2,    0, imax+1, 0, jmax+1);
    free_matrix(_duvdx,     0, imax+1, 0, jmax+1);
    free_matrix(_dv2dy,     0, imax+1, 0, jmax+1);

    /******** CALCULATE G END ********/
    
    for(j=0; j<=jmax+1; j++)
    {
        F[0][j]=U[0][j];
        F[imax+1][j]=U[imax+1][j];
    }
    
    for(i=0; i<=imax+1; i++)
    {
        G[i][0]=V[i][0];
        G[i][jmax+1]=V[i][jmax+1];
    }
}

void calculate_dt(
  double Re,
  double tau,
  double *dt,
  double dx,
  double dy,
  int imax,
  int jmax,
  double **U,
  double **V
)
{
    double umax = 0;
    double vmax = 0;
    int i, j;

    double dtcond, dxcond, dycond;
    double minval;

    /******** Determine umax and vmax *********/
    for(i = 0; i < imax; i++)
    {
        for(j = 0; j < jmax; j++)
        {
            if(umax < U[i][j])
                umax = U[i][j];
            if(vmax < V[i][j])
                vmax = V[i][j];
        }
    }

    /******** Calculate conditions ********/
    dtcond = Re/(2*(1/(dx*dx) + 1/(dy*dy)));
    dxcond = dx/fabs(umax);
    dycond = dy/fabs(vmax);

    /******** Determine smalles condition ********/
    minval = dtcond;
    if(minval > dxcond)
        minval = dxcond;
    if(minval > dycond)
        minval = dycond;
    
    /******** Calculate dt ********/
    *dt = tau * minval;
}

void calculate_uv(
  double dt,
  double dx,
  double dy,
  int imax,
  int jmax,
  double **U,
  double **V,
  double **F,
  double **G,
  double **P
)
{
    int i, j;
    double dtodx, dtody;  
    
    /******** Calculate dt/dx and dt/dy (it is the same for each element) ********/
    dtodx = dt/dx;
    dtody = dt/dy;
    
    /******** Calculate u in step next step ********/
    for(i = 1; i < imax; i++)
    {
        for(j = 1; j < jmax+1; j++)
        {
            U[i][j] = F[i][j] - dtodx*(P[i+1][j] - P[i][j]);
        }
    }

    /******** Calculate v in step next step ********/
    for(i = 1; i < imax + 1; i++)
    {
        for(j = 1; j < jmax; j++)
        {
            V[i][j] = G[i][j] - dtody*(P[i][j+1] - P[i][j]);
        }
    }
}

void calculate_rs(
  double dt,
  double dx,
  double dy,
  int imax,
  int jmax,
  double **F,
  double **G,
  double **RS
)
{
    int i, j;
   
    /******** Calculate RS ********/

    for(i = 1; i < imax + 1; i++)
    {
        for(j = 1; j < jmax + 1; j++)
        {
            RS[i][j]= 1/dt*((F[i][j]-F[i-1][j])/dx+(G[i][j]-G[i][j-1])/dy);
            }
        }
}
