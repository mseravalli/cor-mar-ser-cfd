#include "boundary_val.h"
#include "helper.h"
#include "uvp.h"
#include <math.h>
#include "mpi.h"
#include "parallel.h"

void calculate_fg(
  double Re,
  double GX,
  double GY,
  double alpha,
  double dt,
  double dx,
  double dy,
  int il,
  int ir,
  int jt,
  int jb,
  int rank_l,
  int rank_r,
  int rank_b,
  int rank_t,
  double **U,
  double **V,
  double **F,
  double **G
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

    int istart, iend;
    int jstart, jend;
    
    int imax, jmax;
    imax = ir-il+1;
    jmax = jt-jb+1;
    
    /******** VARIABLE DECLARATION END ********/


    /******** CALCULATE F START ********/
    
    /******** Calculate F ********/

    istart = (rank_l == MPI_PROC_NULL? 2 : 1 );
    iend = (rank_r == MPI_PROC_NULL ? imax : imax + 1);

    for(i = istart; i <= iend; i++)
    {
        for(j = 1; j <= jmax; j++)
        {

            /* du2dx */
            firstOperand = ( 1 / dx) *
                ( pow((U[i][j] + U[i+1][j])/2, 2) -  pow((U[i-1][j] + U[i][j])/2, 2));

            secondOperand = ( alpha / dx ) *
                ( ( (abs(U[i][j] + U[i+1][j])/2) * ((U[i][j] - U[i+1][j])/2) )  
                - ( (abs(U[i-1][j] + U[i][j])/2) * ((U[i-1][j] - U[i][j])/2) ) );

            du2dx = firstOperand + secondOperand;
            
            /* duvdy */
            firstOperand = ( 1 / dy )*
            ( ( (V[i-1][j+1] + V[i][j+1]) / 2 ) * ( (U[i][j] + U[i][j+1]) / 2 )
            - ( (V[i-1][j] + V[i][j]) / 2 ) * ( (U[i][j-1] + U[i][j]) / 2 ) );

            secondOperand = ( alpha / dy )*
            ( ( (abs(V[i-1][j+1] + V[i][j+1])/2) * ((U[i][j] - U[i][j+1])/2) )
            - ( (abs(V[i-1][j] + V[i][j])/2) * ((U[i][j-1] - U[i][j])/2) ));
            duvdy = firstOperand + secondOperand;

            /* d2udx2 */
            d2udx2 = (U[i+1][j] - 2*U[i][j] + U[i-1][j]) / pow(dx,2);

            /* d2udy2 */
            d2udy2 = (U[i][j+1] - 2*U[i][j] + U[i][j-1]) / pow(dy,2);

            F[i][j] = U[i][j] + dt * ( (1/Re) * (d2udx2 + d2udy2 ) - du2dx - duvdy + GX );
        }
    }

    /******** CALCULATE F END ********/

    /******** CALCULATE G START ********/

    /******** Calculate G ********/

    jstart = (rank_b == MPI_PROC_NULL ? 2 : 1 );
    jend = (rank_t == MPI_PROC_NULL ? jmax : jmax + 1);

    for(i = 1; i <= imax; i++)
    {
        for(j = jstart; j <= jend; j++)
        {

           /* dv2dy */ 
            firstOperand = ( 1 / dy) *
                ( pow((V[i][j] + V[i][j+1])/2, 2) -  pow((V[i][j-1] + V[i][j])/2, 2));

            secondOperand = ( alpha / dy ) *
                ( ( (abs(V[i][j] + V[i][j+1])/2) * ((V[i][j] - V[i][j+1])/2) )  
                - ( (abs(V[i][j-1] + V[i][j])/2) * ((V[i][j-1] - V[i][j])/2) ) );

            dv2dy = firstOperand + secondOperand;

            /* duvdx */
            firstOperand = ( 1 / dx )*
            ( ( (U[i+1][j-1] + U[i+1][j]) / 2 ) * ( (V[i][j] + V[i+1][j]) / 2 )
            - ( (U[i][j-1] + U[i][j]) / 2 ) * ( (V[i-1][j] + V[i][j]) / 2 ) );

            secondOperand = ( alpha / dx )*
            ( ( (abs(U[i+1][j-1] + U[i+1][j])/2) * ((V[i][j] - V[i+1][j])/2) )
            - ( (abs(U[i][j-1] + U[i][j])/2) * ((V[i-1][j] - V[i][j])/2) ));

            duvdx = firstOperand + secondOperand;

            /* d2vdx2 */
            d2vdx2 = (V[i+1][j] - 2*V[i][j] + V[i-1][j]) / pow(dx,2);
          
            /* d2vdy2*/
            d2vdy2 = (V[i][j+1] - 2*V[i][j] + V[i][j-1]) / pow(dy,2);

            G[i][j] = V[i][j] + dt * ( (1/Re) * ( d2vdx2 + d2vdy2 ) - duvdx - dv2dy + GY);
        }
    }

    /******** CALCULATE G END ********/
    
    /******** BOUNDARY VALUES START ********/
    for(j=1; j<=jmax; j++)
    {
        F[istart-1][j]=U[istart-1][j];
        F[iend+1][j]=U[iend+1][j];
    }
    
    for(i=1; i<=imax; i++)
    {
        G[i][jstart-1]=V[i][jstart-1];
        G[i][jend+1]=V[i][jend+1];
    }
    /******** BOUNDARY VALUES END ********/
}

void calculate_dt(
  double Re,
  double tau,
  double *dt,
  double dx,
  double dy,
  int il,
  int ir,
  int jt,
  int jb,
  double **U,
  double **V
)
{
    double umax = 0;
    double vmax = 0;
    double globUmax;
    double globVmax;
    int i, j;

    double dtcond, dxcond, dycond;
    double minval;

    int imax, jmax;

    imax = ir-il+1;
    jmax = jt-jb+1;

    /******** Determine umax and vmax *********/
    for(i = 1; i <= imax+1; i++)
    {
        for(j = 1; j <= jmax+1; j++)
        {
            if(umax < fabs(U[i][j]))
                umax = fabs(U[i][j]);
            if(vmax < fabs(V[i][j]))
                vmax = fabs(V[i][j]);
        }
    }
    
    /* determines globan umax and vmax and places it in master thread */
    MPI_Reduce(&umax, &globUmax, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&vmax, &globVmax, 1, MPI_DOUBLE, MPI_MAX, 1, MPI_COMM_WORLD);
    /* broadcasts umax nad vmax */
    MPI_Bcast(&globUmax, 1, MPI_DOUBLE, 2, MPI_COMM_WORLD);
    MPI_Bcast(&globVmax, 1, MPI_DOUBLE, 3, MPI_COMM_WORLD);
    
    
    /******** Calculate conditions ********/
    dtcond = Re/(2*(1/(dx*dx) + 1/(dy*dy)));
    dxcond = dx/fabs(globUmax);
    dycond = dy/fabs(globVmax);

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
  int il,
  int ir,
  int jt, 
  int jb,
  int rank_l,
  int rank_r,
  int rank_b,
  int rank_t,
  double **U,
  double **V,
  double **F,
  double **G,
  double **P
)
{
    int i, j;
    int istart, iend;
    int jstart, jend;
    double dtodx, dtody;  

    int imax, jmax;

    imax = ir-il+1;
    jmax = jt-jb+1;
    
    /******** Calculate dt/dx and dt/dy (it is the same for each element) ********/
    dtodx = dt/dx;
    dtody = dt/dy;
    
    /******** Calculate u in step next step ********/
    istart = (rank_l == MPI_PROC_NULL? 2 : 1 );
    iend = (rank_r == MPI_PROC_NULL ? imax : imax + 1);

    for(i = istart; i <= iend; i++)
    {
        for(j = 1; j <= jmax; j++)
        {
            U[i][j] = F[i][j] - dtodx*(P[i][j] - P[i-1][j]);
        }
    }

    /******** Calculate v in step next step ********/
    jstart = (rank_b == MPI_PROC_NULL ? 2 : 1 );
    jend = (rank_t == MPI_PROC_NULL ? jmax : jmax + 1);

    for(i = 1; i <= imax; i++)
    {
        for(j = jstart; j <= jend; j++)
        {
            V[i][j] = G[i][j] - dtody*(P[i][j] - P[i][j-1]);
        }
    }
}

void calculate_rs(
  double dt,
  double dx,
  double dy,
  int il,
  int ir,
  int jt,
  int jb,
  int rank_l,
  int rank_r,
  int rank_b,
  int rank_t,
  double **F,
  double **G,
  double **RS
)
{
    int i, j, istart, jstart, iend, jend;

    int imax, jmax;

    imax = ir-il+1;
    jmax = jt-jb+1;
   
    /******** Calculate RS ********/
    for(i = 2; i <= imax+1; i++)
    {
        for(j = 2; j <= jmax+1; j++)
        {
            RS[i-1][j-1]= 1/dt*((F[i][j-1]-F[i-1][j-1])/dx+(G[i-1][j]-G[i-1][j-1])/dy);
        }
    }
}
