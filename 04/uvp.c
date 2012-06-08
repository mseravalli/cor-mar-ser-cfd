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

    
    /******** VARIABLE DECLARATION END ********/


    /******** CALCULATE F START ********/
    
    /******** Calculate F ********/

    istart = (rank_l == MPI_PROC_NULL ? 1 : 1);
    iend = (rank_r == MPI_PROC_NULL ? ir-il+1 : ir-il+1);

    for(i = istart; i <= iend; i++)
    {
        for(j = 1; j <= jt-jb+1; j++)
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
            ( ( (V[i-istart+1][j+1] + V[i-istart+1+1][j+1]) / 2 ) * ( (U[i][j] + U[i][j+1]) / 2 )
            - ( (V[i-istart+1][j] + V[i-istart+1+1][j]) / 2 ) * ( (U[i][j-1] + U[i][j]) / 2 ) );

            secondOperand = ( alpha / dy )*
            ( ( (abs(V[i-istart+1][j+1] + V[i-istart+1+1][j+1])/2) * ((U[i][j] - U[i][j+1])/2) )
            - ( (abs(V[i-istart+1][j] + V[i-istart+1+1][j])/2) * ((U[i][j-1] - U[i][j])/2) ));
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

    jstart = (rank_b == MPI_PROC_NULL ? 1 : 1);
    jend = (rank_t == MPI_PROC_NULL ? jt-jb+1 : jt-jb+1);
    
    for(i = 1; i <= ir-il+1; i++)
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
            ( ( (U[i+1][j-jstart+1] + U[i+1][j-jstart+1+1]) / 2 ) * ( (V[i][j] + V[i+1][j]) / 2 )
            - ( (U[i][j-jstart+1] + U[i][j-jstart+1+1]) / 2 ) * ( (V[i-1][j] + V[i][j]) / 2 ) );

            secondOperand = ( alpha / dx )*
            ( ( (abs(U[i+1][j-jstart+1] + U[i+1][j-jstart+1+1])/2) * ((V[i][j] - V[i+1][j])/2) )
            - ( (abs(U[i][j-jstart+1] + U[i][j-jstart+1+1])/2) * ((V[i-1][j] - V[i][j])/2) ));

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
    if(rank_l == MPI_PROC_NULL)
    {  
      for(j=0; j<=jt-jb+2; j++)
      {
        F[0][j]=U[0][j];
      }
    }

    if(rank_r == MPI_PROC_NULL)
    {  
      for(j=0; j<=jt-jb+2; j++)
      {
        F[ir-il+2][j]=U[ir-il+2][j];
      }
    }
    
    
    if(rank_b == MPI_PROC_NULL)
    {
        for(i=0; i<=ir-il+2; i++)
        {
           G[i][0]=V[i][0];
        }
    }
    
    if(rank_t == MPI_PROC_NULL)
    {
      for(i=0; i<=ir-il+2; i++)
      {
          G[i][jt-jb+2]=V[i][jt-jb+2];
      }
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

    /******** Determine umax and vmax *********/
    for(i = 1; i <= ir-il+1; i++)
    {
        for(j = 1; j <= jt-jb+1; j++)
        {
            if(umax < U[i][j])
                umax = U[i][j];
            if(vmax < V[i][j])
                vmax = V[i][j];
        }
    }
    
    /* determines globan umax and vmax and places it in master thread */
    MPI_Reduce(&umax, &globUmax, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&vmax, &globVmax, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    /* broadcasts umax nad vmax */
    MPI_Bcast(&globUmax, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&globVmax, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    
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
    
    /******** Calculate dt/dx and dt/dy (it is the same for each element) ********/
    dtodx = dt/dx;
    dtody = dt/dy;
    
    /******** Calculate u in step next step ********/
    istart = (rank_l == MPI_PROC_NULL ? 1 : 1);
    iend = (rank_r == MPI_PROC_NULL ? ir - il + 2 : ir - il + 2 );
    for(i = istart; i < iend; i++)
    {
        for(j = 1; j < jt-jb+2; j++)
        {
            U[i][j] = F[i][j] - dtodx*(P[i-istart+1][j] - P[i-istart][j]);
        }
    }

    /******** Calculate v in step next step ********/
    jstart = (rank_b == MPI_PROC_NULL ? 1 : 1);
    jend = (rank_t == MPI_PROC_NULL ? jt - jb + 2 : jt - jb +2);
    for(i = 1; i < ir-il+2; i++)
    {
        for(j = jstart; j < jend; j++)
        {
            V[i][j] = G[i][j] - dtody*(P[i][j-jstart+1] - P[i][j-jstart]);
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
   
    /******** Calculate RS ********/
    istart = (rank_l == MPI_PROC_NULL ? 1 : 1);
    iend = (rank_r == MPI_PROC_NULL ? ir-il+1 : ir-il+1);

    jstart = (rank_b == MPI_PROC_NULL ? 1 : 1);
    jend = (rank_t == MPI_PROC_NULL ? jt-jb+1 : jt-jb+1);
 
    for(i = istart; i <= iend; i++)
    {
        for(j = jstart; j <= jend; j++)
        {
            RS[i-istart+1][j-jstart+1]= 1/dt*((F[i+1][j-jstart+1]-F[i][j-jstart+1])/dx+(G[i-istart+1][j+1]-G[i-istart+1][j])/dy);
            }
        }
}
