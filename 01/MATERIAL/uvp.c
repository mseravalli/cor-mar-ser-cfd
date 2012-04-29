#include "uvp.h"
#include "helper.h"
#include "matrix_op.h"
#include "disc.h"
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
    /******** CALCULATE F START ********/

    /******** Calculate the single derivatives ********/
    double** _d2udx2    = d2udx2(U, dx,     0, imax+1, 0, jmax+1);
    double** _d2udy2    = d2udy2(U, dy,     0, imax+1, 0, jmax+1);
    double** _du2dx     = du2dx(U, dx,      0, imax+1, 0, jmax+1);
    double** _duvdy     = duvdy(U, V, dy,   0, imax+1, 0, jmax+1);

    /******** Calculate F ********/
    double** add_d2udx2_d2udy2  = add_mat(_d2udx2, _d2udy2,         0, imax+1, 0, jmax+1);
    double** mult_Re            = mult_scalar(add_d2udx2_d2udy2, Re,0, imax+1, 0, jmax+1);
    double** sub_du2dx          = sub_mat(mult_Re, _du2dx,          0, imax+1, 0, jmax+1);
    double** sub_duvdy          = sub_mat(sub_du2dx, _duvdy,        0, imax+1, 0, jmax+1);
    double** add_GX             = add_scalar(sub_duvdy, GX,         0, imax+1, 0, jmax+1);
    double** mult_dt            = mult_scalar(add_GX, dt,           0, imax+1, 0, jmax+1);
    F                           = add_mat(U, mult_dt,               0, imax+1, 0, jmax+1);

    /******** Free the calulations ********/
    free_matrix(mult_dt,            0, imax+1, 0, jmax+1);
    free_matrix(add_GX,             0, imax+1, 0, jmax+1);
    free_matrix(sub_duvdy,          0, imax+1, 0, jmax+1);
    free_matrix(sub_du2dx,          0, imax+1, 0, jmax+1);
    free_matrix(mult_Re,            0, imax+1, 0, jmax+1);
    free_matrix(add_d2udx2_d2udy2,  0, imax+1, 0, jmax+1);

    /******** Free the derivatives ********/
    free_matrix(_duvdy,     0, imax+1, 0, jmax+1);
    free_matrix(_du2dx,     0, imax+1, 0, jmax+1);
    free_matrix(_d2udy2,    0, imax+1, 0, jmax+1);
    free_matrix(_d2udx2,    0, imax+1, 0, jmax+1);

    /******** CALCULATE F END ********/
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
   
    dtcond = Re/(2*(1/(dx*dx) + 1/(dy*dy)));
    dxcond = dx/fabs(umax);
    dycond = dy/fabs(vmax);
   
    minval = dtcond;
    if(minval > dxcond)
        minval = dxcond;
    if(minval > dycond)
        minval = dycond;
      
    *dt = tau * minval;
}
