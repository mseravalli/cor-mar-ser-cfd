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
    
    /******** VARIABLE DECLARATION START ********/
    
    double** mult_Re;
    double** mult_dt;

    double** _d2udx2;
    double** _d2udy2;
    double** _du2dx ;
    double** _duvdy ;

    double** add_d2udx2_d2udy2;
    double** sub_du2dx        ;
    double** sub_duvdy        ;
    double** add_GX           ;

    double** _d2vdx2;
    double** _d2vdy2;
    double** _duvdx ;
    double** _dv2dy ;

    double** add_d2vdx2_d2vdy2;
    double** sub_duvdx        ;
    double** sub_dv2dy        ;
    double** add_GY           ;

    /******** VARIABLE DECLARATION END ********/


    /******** CALCULATE F START ********/

    /******** Calculate the single derivatives ********/
    _d2udx2 = d2udx2(U, dx,   0, imax+1, 0, jmax+1);
    _d2udy2 = d2udy2(U, dy,   0, imax+1, 0, jmax+1);
    _du2dx  = du2dx(U, dx,    0, imax+1, 0, jmax+1);
    _duvdy  = duvdy(U, V, dy, 0, imax+1, 0, jmax+1);

    /******** Calculate F ********/
    add_d2udx2_d2udy2 = add_mat(_d2udx2, _d2udy2,            0, imax+1, 0, jmax+1);
    mult_Re           = mult_scalar(add_d2udx2_d2udy2, 1/Re, 0, imax+1, 0, jmax+1);
    sub_du2dx         = sub_mat(mult_Re, _du2dx,             0, imax+1, 0, jmax+1);
    sub_duvdy         = sub_mat(sub_du2dx, _duvdy,           0, imax+1, 0, jmax+1);
    add_GX            = add_scalar(sub_duvdy, GX,            0, imax+1, 0, jmax+1);
    mult_dt           = mult_scalar(add_GX, dt,              0, imax+1, 0, jmax+1);
    F                 = add_mat(U, mult_dt,                  0, imax+1, 0, jmax+1);

    /******** Free the calulations ********/
    free_matrix(mult_dt,           0, imax+1, 0, jmax+1);
    free_matrix(add_GX,            0, imax+1, 0, jmax+1);
    free_matrix(sub_duvdy,         0, imax+1, 0, jmax+1);
    free_matrix(sub_du2dx,         0, imax+1, 0, jmax+1);
    free_matrix(mult_Re,           0, imax+1, 0, jmax+1);
    free_matrix(add_d2udx2_d2udy2, 0, imax+1, 0, jmax+1);

    /******** Free the derivatives ********/
    free_matrix(_duvdy,  0, imax+1, 0, jmax+1);
    free_matrix(_du2dx,  0, imax+1, 0, jmax+1);
    free_matrix(_d2udy2, 0, imax+1, 0, jmax+1);
    free_matrix(_d2udx2, 0, imax+1, 0, jmax+1);

    /******** CALCULATE F END ********/
    

    /******** CALCULATE G START ********/

    /******** Calculate the single derivatives ********/
    _d2vdx2 = d2vdx2(U, dx,   0, imax+1, 0, jmax+1);
    _d2vdy2 = d2vdy2(U, dy,   0, imax+1, 0, jmax+1);
    _duvdx  = duvdx(U, V, dy, 0, imax+1, 0, jmax+1);
    _dv2dy  = dv2dy(U, dx,    0, imax+1, 0, jmax+1);

    /******** Calculate G ********/
    add_d2vdx2_d2vdy2 = add_mat(_d2vdx2, _d2vdy2,            0, imax+1, 0, jmax+1);
    mult_Re           = mult_scalar(add_d2vdx2_d2vdy2, 1/Re, 0, imax+1, 0, jmax+1);
    sub_duvdx         = sub_mat(mult_Re, _duvdx,             0, imax+1, 0, jmax+1);
    sub_dv2dy         = sub_mat(sub_duvdx, _dv2dy,           0, imax+1, 0, jmax+1);
    add_GY            = add_scalar(sub_dv2dy, GY,            0, imax+1, 0, jmax+1);
    mult_dt           = mult_scalar(add_GY, dt,              0, imax+1, 0, jmax+1);
    G                 = add_mat(V, mult_dt,                  0, imax+1, 0, jmax+1);

    /******** Free the calulations ********/
    free_matrix(mult_dt,            0, imax+1, 0, jmax+1);
    free_matrix(add_GY,             0, imax+1, 0, jmax+1);
    free_matrix(sub_dv2dy,          0, imax+1, 0, jmax+1);
    free_matrix(sub_duvdx,          0, imax+1, 0, jmax+1);
    free_matrix(mult_Re,            0, imax+1, 0, jmax+1);
    free_matrix(add_d2vdx2_d2vdy2,  0, imax+1, 0, jmax+1);

    /******** Free the derivatives ********/
    free_matrix(_dv2dy,     0, imax+1, 0, jmax+1);
    free_matrix(_duvdx,     0, imax+1, 0, jmax+1);
    free_matrix(_d2vdy2,    0, imax+1, 0, jmax+1);
    free_matrix(_d2vdx2,    0, imax+1, 0, jmax+1);

    /******** CALCULATE G END ********/
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
