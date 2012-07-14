#include "boundary_val.h"
#include "helper.h"
#include "uvp.h"
#include <math.h>
#include "constants.h"
#include <omp.h>

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
  double **T,
  double beta 
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

    #pragma omp parallel private (i, j, d2udx2, d2udy2, du2dx, duvdy, d2vdx2, d2vdy2, duvdx, dv2dy, firstOperand, secondOperand)
    {
        /******** VARIABLE DECLARATION END ********/


        /******** CALCULATE F START ********/
        
        /******** Calculate F ********/
        
        #pragma omp for collapse(2)
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

                    F[i][j] = U[i][j] + dt * ( (1/Re) * (d2udx2 + d2udy2 ) - du2dx - duvdy + GX*(1-beta*(T[i][j]+T[i+1][j])/2) );
                }
                
            }
        }
        
        /******** CALCULATE F END ********/
        

        /******** CALCULATE G START ********/

        /******** Calculate G ********/

        #pragma omp for collapse(2)
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

                    G[i][j] = V[i][j] + dt * ( (1/Re) * ( d2vdx2 + d2vdy2 ) - duvdx - dv2dy + GY*(1-beta*(T[i][j]+T[i][j+1])/2));
                }
            }
        }

        /******** CALCULATE G END ********/
        
        /******** BOUNDARY VALUES START ********/
        #pragma omp for
        for(j=1; j<=jmax; j++)
        {
            F[0][j]=U[0][j];
            F[imax][j]=U[imax][j];
        }
        
        #pragma omp for
        for(i=1; i<=imax; i++)
        {
            G[i][0]=V[i][0];
            G[i][jmax]=V[i][jmax];
        }
        /******** BOUNDARY VALUES END ********/

    }
}

void calculate_dt(
  double Re,
  double Pr,
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

    double dttcond, dtdcond, dtcond, dxcond, dycond;
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
    dttcond = Re*Pr/(2*(1/(dx*dx) + 1/(dy*dy)));
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

    #pragma omp parallel for private(i, j, k, ducdx, d2cdx2, dvcdy, d2cdy2, firstOperand, secondOperand) collapse(3)
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
  int **Flag,
  double Ei,
  double Er,
  double **T,
  int reactantsNum,
  int productsNum,
  double catRate
)
{
    int i;
    int j;
    int k;
    int l;
    int cat;
    double forwardReaction;
    double backwardReaction;
    double foFactor;
    double bwFactor;

    /****** Calculate Q start ******/
    
    #pragma omp parallel for private(i, j, k, l, cat, forwardReaction, backwardReaction, foFactor, bwFactor) collapse(3)
    for (i = 1; i <= imax; i++) {
        for (j = 1; j <= jmax; j++) {
            for (k = 0; k < kmax; k++){
                if (Flag[i][j] >= C_F) 
		{
                    /* to avoid an additional if:
                     * see if it a catalyst, if it is the result of & will be 32
                     * shift 8 positions to obtain 1 or 0
                     * multiply by the specified catalyst rate
                     */
		    cat = Flag[i][j] & 0x20;
		    cat = cat >> 5;
		    cat = cat*(catRate-1) + 1;
		    foFactor = 1;
		    for(l = 0; l < reactantsNum; l++)
		    {
		        foFactor *= pow(C[l][i][j], -K[0][l]);
	            }
                    forwardReaction = (cat*K[1][k]*exp(-Ei/T[i][j])*foFactor);
		    bwFactor = 1;
		    for(l = reactantsNum; l < reactantsNum + productsNum; l++)
		    {
		    	bwFactor *= pow(C[l][i][j], K[0][l]);
		    }
                    backwardReaction = (K[2][k]*exp(-Er/T[i][j])*bwFactor);
                    Q[k][i][j] = forwardReaction -backwardReaction;
                }
            
            }
        }
    }

}

void calculate_t(
  double dt,
  double dx,
  double dy,
  double alpha,
  int imax,
  int jmax,
  double Re,
  double Pr,
  double **U,
  double **V,
  double **T,
  int **Flag
){
    double dutdx;
    double d2tdx2;
    double dvtdy;
    double d2tdy2;
    double firstOperand;
    double secondOperand;

    int i;
    int j;
    int k;

        for (i = 1; i <= imax; i++) {
            for (j = 1; j <= jmax; j++) {

                if (Flag[i][j] >= C_F) {

                    /*** dutdx ***/

                    firstOperand = (1/dx)*(U[i][j]*(T[i][j]+T[i+1][j])/2 - U[i-1][j]*(T[i-1][j]+T[i][j])/2);
                    secondOperand = (alpha/dx)*(fabs(U[i][j])*(T[i][j]-T[i+1][j])/2 - fabs(U[i-1][j])*(T[i-1][j]-T[i][j])/2);

                    dutdx = firstOperand + secondOperand;

                    /*** dvtdy ***/

                    firstOperand = (1/dy)*(V[i][j]*(T[i][j]+T[i][j+1])/2 - V[i][j-1]*(T[i][j-1]+T[i][j])/2);
                    secondOperand = (alpha/dy)*(fabs(V[i][j])*(T[i][j]-T[i][j+1])/2 - fabs(V[i][j-1])*(T[i][j-1]-T[i][j])/2);

                    dvtdy = firstOperand + secondOperand;

                    /*** d2tdx2 ***/

                    d2tdx2 = (T[i+1][j]-2*T[i][j]+T[i-1][j])/(dx*dx);

                    /*** d2tdy2 ***/

                    d2tdy2 = (T[i][j+1]-2*T[i][j]+T[i][j-1])/(dy*dy);

                    /*** T(t+dt) ***/

                    T[i][j] =  T[i][j] + dt*(1/(Re*Pr)*(d2tdx2+d2tdy2) - dutdx - dvtdy);
                }
            }
        }

}
