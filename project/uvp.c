#include "boundary_val.h"
#include "helper.h"
#include "uvp.h"
#include <math.h>
#include "constants.h"

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
  double **G,
  int **Flag,
  int wl,
  int wr,
  int wt,
  int wb
) 
{
    
    /******** VARIABLE DECLARATION START ********/

    double d2udx2;
    double d2udy2;
    double du2dx ;
    double duvdy ;

    double d2vdx2;
    double d2vdy2;
    double duvdx ;
    double dv2dy ;
    
    double firstOperand;
    double secondOperand;

    int i;
    int j;

    /******** VARIABLE DECLARATION END ********/


    /******** CALCULATE F START ********/
    
    /******** Calculate F ********/
    
    for(i = 1; i <= imax-1; i++)
    {
        for(j = 1; j <= jmax; j++)
        {

            /* Fluid */
            if (Flag[i][j]>= C_F && Flag[i+1][j]>= C_F)
            {
                /* du2dx */
                firstOperand = ( 1 / dx) *
                    ( ((U[i][j] + U[i+1][j])/2)*((U[i][j] + U[i+1][j])/2) - ((U[i-1][j] + U[i][j])/2)*((U[i-1][j] + U[i][j])/2));

                secondOperand = ( alpha / dx ) *
                    ( ( (abs(U[i][j] + U[i+1][j])/2) * ((U[i][j] - U[i+1][j])/2) )  
                    - ( (abs(U[i-1][j] + U[i][j])/2) * ((U[i-1][j] - U[i][j])/2) ) );

                du2dx = firstOperand + secondOperand;
                
                /* duvdy */
                firstOperand = ( 1 / dy )*
                ( ( (V[i][j] + V[i+1][j]) / 2 ) * ( (U[i][j] + U[i][j+1]) / 2 )
                - ( (V[i][j-1] + V[i+1][j-1]) / 2 ) * ( (U[i][j-1] + U[i][j]) / 2 ) );

                secondOperand = ( alpha / dy )*
                ( ( (abs(V[i][j] + V[i+1][j])/2) * ((U[i][j] - U[i][j+1])/2) )
                - ( (abs(V[i][j-1] + V[i+1][j-1])/2) * ((U[i][j-1] - U[i][j])/2) ));
                duvdy = firstOperand + secondOperand;

                /* d2udx2 */
                d2udx2 = (U[i+1][j] - 2*U[i][j] + U[i-1][j]) / (dx*dx);

                /* d2udy2 */
                d2udy2 = (U[i][j+1] - 2*U[i][j] + U[i][j-1]) / (dy*dy);

                F[i][j] = U[i][j] + dt * ( (1/Re) * (d2udx2 + d2udy2 ) - du2dx - duvdy + GX );
            }
            
        }
    }
    
    /******** CALCULATE F END ********/
    

    /******** CALCULATE G START ********/

    /******** Calculate G ********/

    for(i = 1; i <= imax; i++)
    {
        for(j = 1; j <= jmax-1; j++)
        {
            /* Fluid */
            if (Flag[i][j]>=C_F && Flag[i][j+1]>=C_F)
            {
               /* dv2dy */ 
                firstOperand = ( 1 / dy) *
                    ( ((V[i][j] + V[i][j+1])/2)*((V[i][j] + V[i][j+1])/2) - ((V[i][j-1] + V[i][j])/2)*((V[i][j-1] + V[i][j])/2));

                secondOperand = ( alpha / dy ) *
                    ( ( (abs(V[i][j] + V[i][j+1])/2) * ((V[i][j] - V[i][j+1])/2) )  
                    - ( (abs(V[i][j-1] + V[i][j])/2) * ((V[i][j-1] - V[i][j])/2) ) );

                dv2dy = firstOperand + secondOperand;

                /* duvdx */
                firstOperand = ( 1 / dx )*
                ( ( (U[i][j] + U[i][j+1]) / 2 ) * ( (V[i][j] + V[i+1][j]) / 2 )
                - ( (U[i-1][j] + U[i-1][j+1]) / 2 ) * ( (V[i-1][j] + V[i][j]) / 2 ) );

                secondOperand = ( alpha / dx )*
                ( ( (abs(U[i][j] + U[i][j+1])/2) * ((V[i][j] - V[i+1][j])/2) )
                - ( (abs(U[i-1][j] + U[i-1][j+1])/2) * ((V[i-1][j] - V[i][j])/2) ));

                duvdx = firstOperand + secondOperand;

                /* d2vdx2 */
                d2vdx2 = (V[i+1][j] - 2*V[i][j] + V[i-1][j]) / (dx*dx);
              
                /* d2vdy2*/
                d2vdy2 = (V[i][j+1] - 2*V[i][j] + V[i][j-1]) / (dy*dy);

                G[i][j] = V[i][j] + dt * ( (1/Re) * ( d2vdx2 + d2vdy2 ) - duvdx - dv2dy + GY);
            }
        }
    }

    /******** CALCULATE G END ********/
    
    /******** BOUNDARY VALUES START ********/
    for(j=1; j<=jmax; j++)
    {
        F[0][j]=U[0][j];
        F[imax][j]=U[imax][j];
    }
    
    for(i=1; i<=imax; i++)
    {
        G[i][0]=V[i][0];
        G[i][jmax]=V[i][jmax];
    }
    /******** BOUNDARY VALUES END ********/
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
  double **V,
  int **Flag,
  double D
)
{
    double umax = 0;
    double vmax = 0;
    int i, j;

    double dtdcond, dtcond, dxcond, dycond;
    double minval;

    /******** Determine umax and vmax *********/
    for(i = 0; i < imax; i++)
    {
        for(j = 0; j < jmax; j++)
        {
            if(Flag[i][j] >= C_F){
                if(umax < U[i][j])
                    umax = U[i][j];
                if(vmax < V[i][j])
                    vmax = V[i][j];
            }
        }
    }

    /******** Calculate conditions ********/
    dtdcond = D/(2*(1/(dx*dx) + 1/(dy*dy)));
    dtcond = Re/(2*(1/(dx*dx) + 1/(dy*dy)));
    dxcond = dx/fabs(umax);
    dycond = dy/fabs(vmax);

    /******** Determine smalles condition ********/
    minval = dtdcond;
    if(minval > dtcond)
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
  double **P,
  int **Flag
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
            if (Flag[i][j] >= C_F && Flag[i+1][j] >= C_F){
                 U[i][j] = F[i][j] - dtodx*(P[i+1][j] - P[i][j]);
            }
        }
    }

    /******** Calculate v in step next step ********/
    for(i = 1; i < imax + 1; i++)
    {
        for(j = 1; j < jmax; j++)
        {
            if (Flag[i][j] >= C_F && Flag[i][j+1] >= C_F){
                V[i][j] = G[i][j] - dtody*(P[i][j+1] - P[i][j]);
            }
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
  double **RS,
  int **Flag
)
{
    int i, j;
   
    /******** Calculate RS ********/

    for(i = 1; i <= imax; i++)
    {
        for(j = 1; j <= jmax; j++)
        {
            if (Flag[i][j] >= C_F){
                RS[i][j]= 1/dt*((F[i][j]-F[i-1][j])/dx+(G[i][j]-G[i][j-1])/dy);
            }
        }
     }
}

void calculate_c(
  double dt,
  double dx,
  double dy,
  double alpha,
  double D,
  int imax,
  int jmax,
  int kmax,
  double **U,
  double **V,
  double ***Q,
  double ***C,
  int **Flag
)
{
    double ducdx;
    double d2cdx2;
    double dvcdy;
    double d2cdy2;
    double firstOperand;
    double secondOperand;

    int i;
    int j;
    int k;

    /****** Calculate C start ******/

    for (k = 0; k < kmax; k++){
        for (i = 1; i <= imax; i++) {
            for (j = 1; j <= jmax; j++) {

                if (Flag[i][j] >= C_F) {

                    /*** ducdx ***/

                    firstOperand = (1/dx)*(U[i][j]*(C[k][i][j]+C[k][i+1][j])/2 - U[i-1][j]*(C[k][i-1][j]+C[k][i][j])/2);
                    secondOperand = (alpha/dx)*(fabs(U[i][j])*(C[k][i][j]-C[k][i+1][j])/2 - fabs(U[i-1][j])*(C[k][i-1][j]-C[k][i][j])/2);

                    ducdx = firstOperand + secondOperand;

                    /*** dvcdy ***/

                    firstOperand = (1/dy)*(V[i][j]*(C[k][i][j]+C[k][i][j+1])/2 - V[i][j-1]*(C[k][i][j-1]+C[k][i][j])/2);
                    secondOperand = (alpha/dy)*(fabs(V[i][j])*(C[k][i][j]-C[k][i][j+1])/2 - fabs(V[i][j-1])*(C[k][i][j-1]-C[k][i][j])/2);

                    dvcdy = firstOperand + secondOperand;

                    /*** d2cdx2 ***/

                    d2cdx2 = (C[k][i+1][j]-2*C[k][i][j]+C[k][i-1][j])/(dx*dx);

                    /*** d2cdy2 ***/

                    d2cdy2 = (C[k][i][j+1]-2*C[k][i][j]+C[k][i][j-1])/(dy*dy);

                    /*** C(t+dt) ***/

                    C[k][i][j] =  C[k][i][j] + dt*(D*(d2cdx2+d2cdy2) - ducdx - dvcdy + Q[k][i][j]);
                }
            }
        }
    }
}

void calculate_q(
  double **K,
  int imax,
  int jmax,
  int kmax,
  double ***Q,
  double ***C,
  int **Flag
)
{
    int i;
    int j;
    int k;


    /****** Calculate Q start ******/

    /* TODO: build this function again !!*/

    for (i = 1; i <= imax; i++) {
        for (j = 1; j <= jmax; j++) {
            for (k = 0; k < kmax; k++){
                if (Flag[i][j] >= C_F) {
                    Q[k][i][j] = K[1][k]*pow(C[0][i][j],-K[0][0])*pow(C[1][i][j],-K[0][1]);
                   /* Q[1][i][j] = -K[1]*pow(C[0][i][j],K[0])*pow(C[1][i][j],K[1]);
                    Q[2][i][j] =  K[2]*pow(C[0][i][j],K[0])*pow(C[1][i][j],K[1]);
                    Q[3][i][j] =  K[3]*pow(C[0][i][j],K[0])*pow(C[1][i][j],K[1]);*/
                }
            
            }
        }
    }

}
