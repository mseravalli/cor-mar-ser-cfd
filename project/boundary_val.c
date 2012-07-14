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
  int wb,
  int wt,
  int cl,
  int cr,
  int cb,
  int ct,
  double tl,
  double tr,
  double tb,
  double tt,
  double** U,
  double** V,
  double** F,
  double** G,
  double** P,
  double*** C,
  double** T,
  int** Flag,
  int** Sources,
  double* C0
)
{
    int i = 0;
    int j = 0;
    int k = 0;

    #pragma omp parallel private(i, j, k)
    {
        /******** EXTERNAL WALLS START ********/

        /**** VELOCITIES START ****/
        #pragma omp sections
        {
            /** left wall **/
            #pragma omp section
            {
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
            }

            /** right wall **/
            #pragma omp section
            {
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
            }

            /** lower wall **/
            #pragma omp section
            {
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

            /** upper wall **/
            #pragma omp section
            {
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
            }
        }
        /**** VELOCITIES END ****/

        /**** CONCENTRATION START ****/

        /* 
         * if cx specifies a release of a concentration use it 
         * otherwise use the default behaviour for that boundary i.e
         * left boundary no entering concentration => dirichlet with 0
         * right boundary neumann 0
         * bottom boundary neumann 0
         * top boundary neumann 0
         */
        for (k = 0; k < kmax; ++k) {

            /** left and right wall **/
            #pragma omp for
            for (j = 0; j <= jmax; ++j) {
                if (cl & (1 << k)) {
                    C[k][0][j] = 2*C0[k]-C[k][1][j];
                } else {
                    C[k][0][j] = -C[k][1][j];
                }

                if (cr & (1<<k)) {
                    C[k][imax+1][j] = 2*C0[k]-C[k][imax][j];
                } else {
                    C[k][imax+1][j] = C[k][imax][j];
                }
            }

            /** lower and upper wall **/
            #pragma omp for 
            for(i = 0; i <= imax; ++i) {
                if (cb & (1<<k)) {
                    C[k][i][0] = 2*C0[k] - C[k][i][1];
                } else {
                    C[k][i][0] = C[k][i][1];
                }

                if (ct & (1<<k)) {
                    C[k][i][jmax+1] = 2*C0[k] - C[k][i][jmax];
                } else {
                    C[k][i][jmax+1] = C[k][i][jmax];
                }
            }
        }

        /**** CONCENTREATION END ****/

	/**** TEMPERATURE START ****/
        /* if the temperature at a wall is 0 use neumann boundary conditions
         * for neumann bc we use only the case that the bounday has a perfect
         * absorbtion of the heat ie dT/dn = 0
         * otherwise use dirichlet with the given temperature
         */
	
	/** left and right wall **/
	#pragma omp for
	for ( j = 1; j <= jmax; ++j ) {
            /* left wall */
            if (tl > 0) {
                T[0][j] = 2*tl - T[1][j];
            } else {
		double dTdn = 0;
                T[0][j] = T[1][j] + (dTdn * (j - 0.5) * dy * dx);
            }
            /* right wall */
            if (tr > 0) {
                T[imax+1][j] = 2*tr - T[imax][j];
            } else {
		double dTdn = 0;
                T[imax+1][j] = T[imax][j] + (dTdn * (j - 0.5) * dy * dx);
            }
	}

	/** bottom and top wall **/
	#pragma omp for
	for (i = 1; i <= imax; ++i) {
            /* bottom wall */
            if (tb > 0) {
		T[i][0] = 2*tb - T[i][1];
            } else {
		double dTdn = 0;
                T[i][0] = T[i][1] + (dTdn*(i - 0.5) *dx * dy);
            }

            /* top wall */
            if (tt > 0) {
		T[i][jmax+1] = 2*tt - T[i][jmax];
            } else {
		double dTdn = 0;
                T[i][jmax+1] = T[i][jmax] + (dTdn * (i-0.5)*dx*dy);
            }
	}

	/**** TEMPERATURE END ****/

        
        /******** EXTERNAL WALLS END ********/


        /******** INTERNAL OBJECTS START ********/

        /**** VELOCITIES START ****/

        #pragma omp for collapse(2)
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

        /**** VELOCITIES END ****/

        /**** CONCENTRATIONS START ****/
        
        #pragma omp for collapse(3)
        for (k = 0; k < kmax; ++k) {
            for (i = 1; i <= imax; ++i) { 
                for (j = 1; j <= jmax; ++j) {
                    
                    /** if source then dirichlet boundary conditions **/
                    /* shift 1 k times to see if the obstacle is source for
                     * that concentration
                     */
                    if( Sources[i][j] & (1 << k)) {
                        if (Flag[i][j] == B_N)
                        {
                            C[k][i][j] = 2*C0[k]-C[k][i][j+1];
                        }
                        
                        else if (Flag[i][j] == B_W)
                        {
                            C[k][i][j] = 2*C0[k]-C[k][i-1][j];
                        }

                        else if (Flag[i][j] == B_S)
                        {
                            C[k][i][j] = 2*C0[k]-C[k][i][j-1];
                        }

                        else if (Flag[i][j] == B_O)
                        {
                            C[k][i][j] = 2*C0[k]-C[k][i+1][j];
                        }

                        else if (Flag[i][j] == B_NO) 
                        {
                            C[k][i][j] = 2*C0[k]-(C[k][i+1][j]+C[k][i][j+1])/2;
                        }

                        else if (Flag[i][j] == B_NW) 
                        {
                            C[k][i][j] = 2*C0[k]-(C[k][i-1][j]+C[k][i][j+1])/2;
                        }

                        else if (Flag[i][j] == B_SO) 
                        {
                            C[k][i][j] = 2*C0[k]-(C[k][i+1][j]+C[k][i][j-1])/2;
                        }

                        else if (Flag[i][j] == B_SW) 
                        {
                            C[k][i][j] = 2*C0[k]-(C[k][i-1][j]+C[k][i][j-1])/2;
                        }
                    }
                    /** if not source then neumann boundary conditions **/
                    else {
                        if (Flag[i][j] == B_N)
                        {
                            C[k][i][j]=C[k][i][j+1];
                        }
                        
                        else if (Flag[i][j] == B_W)
                        {
                            C[k][i][j]=C[k][i-1][j];
                        }

                        else if (Flag[i][j] == B_S)
                        {
                            C[k][i][j]=C[k][i][j-1];
                        }

                        else if (Flag[i][j] == B_O)
                        {
                            C[k][i][j]=C[k][i+1][j];
                        }

                        else if (Flag[i][j] == B_NO) 
                        {
                            C[k][i][j]=(C[k][i+1][j]+C[k][i][j+1])/2;
                        }

                        else if (Flag[i][j] == B_NW) 
                        {
                            C[k][i][j]=(C[k][i-1][j]+C[k][i][j+1])/2;
                        }

                        else if (Flag[i][j] == B_SO) 
                        {
                            C[k][i][j]=(C[k][i+1][j]+C[k][i][j-1])/2;
                        }

                        else if (Flag[i][j] == B_SW) 
                        {
                            C[k][i][j]=(C[k][i-1][j]+C[k][i][j-1])/2;
                        }
                    }
                }
            }
        }

        /**** CONCENTRATIONS END ****/

	/**** TEMPERATURE START ****/
        #pragma omp for collapse(2)
	for (i = 1; i <= imax; ++i) {
            for (j = 1; j <= jmax; ++j) {
                if (Flag[i][j] == B_N)
                {
                    T[i][j]=T[i][j+1];
                }
                
                else if (Flag[i][j] == B_W)
                {
                    T[i][j]=T[i-1][j];
                }

                else if (Flag[i][j] == B_S)
                {
                    T[i][j]=T[i][j-1];
                }

                else if (Flag[i][j] == B_O)
                {
                    T[i][j]=T[i+1][j];
                }

                else if (Flag[i][j] == B_NO) 
                {
                    T[i][j]=(T[i+1][j]+T[i][j+1])/2;
                }

                else if (Flag[i][j] == B_NW) 
                {
                    T[i][j]=(T[i-1][j]+T[i][j+1])/2;
                }

                else if (Flag[i][j] == B_SO) 
                {
                    T[i][j]=(T[i+1][j]+T[i][j-1])/2;
                }

                else if (Flag[i][j] == B_SW) 
                {
                    T[i][j]=(T[i-1][j]+T[i][j-1])/2;
                }
            }
        }

	/**** TEMPERATURE END ****/

        /******** INTERNAL OBJECTS END ********/
    }
 
}

void spec_boundary_val(
    char *problem,
    int imax,
    int jmax,
    int kmax,
    double dx,
    double dy,
    double Re,
    double deltaP,
    double **U,
    double **V,
    double **P,
    double*** C
    )
{
    int i;
    int j;
    int k;
    /*
    double C0 = 0;
    double C1 = 0;
    */
    /*double dpdx;*/

    if (strcmp(problem, "karman") == 0) {
        for(j = 1; j <= jmax; ++j)
        {
            U[0][j] = 1.0;
            V[0][j] = -V[1][j];
        }
    }

    if (strcmp(problem, "diffusion") == 0) {
        for(j = 1; j <= jmax; ++j)
        {
            U[0][j] = 0;
            V[0][j] = 0;
	    for(k=0; k < kmax; k++)
	    {
	        /* neumann on left wall */
                C[k][0][j] = C[k][1][j];
	    }
        }
    }

    if (strcmp(problem, "semibaffle")==0) {
        for(i = 1; i <= imax; ++i ) {
            U[i][0] = -2 - U[i][1];
            U[i][jmax+1] = -2 - U[i][jmax];
        }
    }

/*
    for (i = 1; i <= imax; ++i) {
        C[0][i][0] = 2*C0-C[0][i][1];
    }

    for (j = 1; j <= jmax; ++j) {
        C[1][0][j] = 2*C1-C[1][1][j];
    }
    */

}

