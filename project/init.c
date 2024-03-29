#include "helper.h"
#include "init.h"
#include "constants.h"

int read_parameters( const char *szFileName,       /* name of the file */
                    double *Re,                /* reynolds number   */
                    double *Pr,
                    double *UI,                /* velocity x-direction */
                    double *VI,                /* velocity y-direction */
                    double *PI,                /* pressure */
                    double *TI,
                    double *GX,                /* gravitation x-direction */
                    double *GY,                /* gravitation y-direction */
                    double *t_end,             /* end time */
                    double *xlength,           /* length of the domain x-dir.*/
                    double *ylength,           /* length of the domain y-dir.*/
                    double *dt,                /* time step */
                    double *dx,                /* length of a cell x-dir. */
                    double *dy,                /* length of a cell y-dir. */
                    int  *imax,                /* number of cells x-direction*/
                    int  *jmax,                /* number of cells y-direction*/
                    double *alpha,             /* uppwind differencing factor*/
                    double *beta,
                    double *omg,               /* relaxation factor */
                    double *tau,               /* safety factor for time step*/
                    int  *itermax,             /* max. number of iterations  */
		                                       /* for pressure per time step */
                    double *eps,               /* accuracy bound for pressure*/
                    int    *wl,
                    int    *wr,
                    int    *wb,
                    int    *wt,
                    int    *cl,
                    int    *cr,
                    int    *cb,
                    int    *ct,
                    double *tl,
                    double *tr,
                    double *tb,
                    double *tt,
                    double *dt_value,           /* time for output */
                    double *deltaP,
                    double* D,
                    int    *kmax,
                    double *ki,                 /* Kinetic of the irreversible reaction */
                    double *kr,                 /* Kinetic of the reversible reaction */
                    double *Ei,
                    double *Er,
                    double* catRate
) 
{
   READ_DOUBLE( szFileName, *xlength );
   READ_DOUBLE( szFileName, *ylength );

   READ_DOUBLE( szFileName, *Re    );
   READ_DOUBLE( szFileName, *Pr    );
   READ_DOUBLE( szFileName, *t_end );
   READ_DOUBLE( szFileName, *dt    );

/* those should be now read from the image
   READ_INT   ( szFileName, *imax );
   READ_INT   ( szFileName, *jmax );
   */

   READ_DOUBLE( szFileName, *omg   );
   READ_DOUBLE( szFileName, *eps   );
   READ_DOUBLE( szFileName, *tau   );
   READ_DOUBLE( szFileName, *alpha );
   READ_DOUBLE( szFileName, *beta  );
   
   READ_INT( szFileName, *wl );
   READ_INT( szFileName, *wr );
   READ_INT( szFileName, *wt );
   READ_INT( szFileName, *wb );

   READ_INT( szFileName, *cl );
   READ_INT( szFileName, *cr );
   READ_INT( szFileName, *ct );
   READ_INT( szFileName, *cb );

   READ_DOUBLE( szFileName, *tl );
   READ_DOUBLE( szFileName, *tr );
   READ_DOUBLE( szFileName, *tt );
   READ_DOUBLE( szFileName, *tb );

   READ_INT   ( szFileName, *itermax );
   READ_DOUBLE( szFileName, *dt_value );

   READ_DOUBLE( szFileName, *UI );
   READ_DOUBLE( szFileName, *VI );
   READ_DOUBLE( szFileName, *GX );
   READ_DOUBLE( szFileName, *GY );
   READ_DOUBLE( szFileName, *PI );
   READ_DOUBLE( szFileName, *TI );

   READ_DOUBLE( szFileName, *deltaP );
   READ_DOUBLE( szFileName, *D );

   READ_INT( szFileName, *kmax );
   READ_DOUBLE( szFileName, *ki );
   READ_DOUBLE( szFileName, *kr );
   READ_DOUBLE( szFileName, *Ei );
   READ_DOUBLE( szFileName, *Er );
   READ_DOUBLE( szFileName, *catRate );



/*
   *dx = *xlength / (double)(*imax);
   *dy = *ylength / (double)(*jmax);
*/
   return 1;
}

void init_C0K(const char *szFileName, int kmax, double* C0, double** K, double ki, double kr, int *reactantsNum, int *productsNum) {
    
    int k;
    char* baseNameC0 = "C";
    char* baseNameK = "K";
    char varNameC0[64];
    char varNameK[64];
    double var;

    for (k = 0; k < kmax; ++k){
        sprintf(varNameC0, "%s%d", baseNameC0, k);
        sprintf(varNameK, "%s%d", baseNameK, k);
        read_double(szFileName, varNameC0, &var);
        C0[k] = var; 
        read_double(szFileName, varNameK, &var);
        K[0][k] = var; 
    }

    *reactantsNum = 0;
    *productsNum = 0;

    for(k = 0; k < kmax; ++k)
    {
    	if(K[0][k] < 0)
	{
	    (*reactantsNum)++;
	}
	else
	{
	    (*productsNum)++;
	}
    }

    for (k = 0; k < kmax; k++){
        K[1][k] = -ki*K[0][k]/K[0][0];
        K[2][k] = kr*K[0][k]/K[0][*reactantsNum];
    }

}

/**
 * The arrays U,V and P are initialized to the constant values UI, VI and PI on
 * the whole domain.
 */
void init_uvp(double UI,
              double VI,
              double PI,
              double TI,
              int imax,
              int jmax,
              char* problem,
              double **U,
              double **V,
              double **P,
              double*** C,
              double **T,
              int kmax)
{
    
    int i ,j;
    int k;

    /* initialize the matrices */
    init_matrix(U, 0, imax + 1, 0, jmax + 1, UI);
    if (strcmp(problem, "step") == 0) {
        for (i = 0; i <=imax + 1; ++i) {
            for (j = 0; j < jmax/2; ++j) {
                U[i][j] = 0;
            }
        }
    }
    
    init_matrix(V, 0, imax + 1, 0, jmax + 1, VI);
    init_matrix(T, 1, imax, 1, jmax, TI);
    init_matrix(P, 1, imax, 1, jmax, PI);

    for (k = 0; k < kmax; k++){
        init_matrix(C[k], 1, imax, 1, jmax, 0);
    }

}

/**
 * Initializes flag matrix
**/
int init_flag(
    int **Problem,
    int imax,
    int jmax,
    int **Flag,
    int **Sources)
{
    int i, j;

    /* 
     * 255 - fluid
     * 128 - catalyst
     * 0 - obstacle
     * 1=(01)b - source of the first substance, 2=(10)b - source of the secont substance, 3=(11)b - source of both substances
     */
    

    for (i = 1; i < imax+1; i++) {
        for (j = 1; j < jmax+1; j++) {
            if(Problem[i][j] != 0 && Problem[i][j] != 255 && Problem[i][j] != 128)
            {
                Sources[i][j] = Problem[i][j];
            }
            else
            {
                Sources[i][j] = 0;
            }
	    }
    }

    for (i = 1; i < imax+1; i++) {
        for (j = 1; j < jmax+1; j++) {
            /*if fluid Problem is 1, Flag set*/
            if(Problem[i][j] == 255)
            {
                Problem[i][j] = 1;
		Flag[i][j] = C_F;
            }
	    /*if catalyst Problem is also 1, Flag set*/
	    else if(Problem[i][j] == 128)
	    {
	    	Problem[i][j] = 1;
		Flag[i][j] = C_C;
	    }
	    /*if obstacle Problem is 0*/
            else
            {
                Problem[i][j] = 0;
            }
        }
    }
    
    for(i = 1; i < imax+1; i++)
    {
        for(j = 1; j < jmax+1; j++)
        {
	    /*Flag for fluid and catalyst is set in the previous loop, so here we are dealing only with obstacles*/
	    if(Problem[i][j] == 0)
            {
                Flag[i][j] = 8 * Problem[i+1][j] + 4 * Problem[i-1][j] + 2 * Problem[i][j-1] + 1 * Problem[i][j+1];
                /*if falg is not valid it returns a wrong result*/
                if(Flag[i][j] == 3 || Flag[i][j] == 7 || Flag[i][j] == 11 || Flag[i][j] == 12 || Flag[i][j] == 13 || Flag[i][j] == 14 || Flag[i][j] == 15){
                    return 1;
                }
            }
        }
    }
    
    /*boundaries depend on only one neighbouring cell*/
    for(i = 1; i < imax+1; i++)
    {
        Flag[i][0] = 8 * Problem[i][1];
        Flag[i][jmax+1] = 4 * Problem[i][jmax];
    }
    
    for(j = 1; j < jmax+1; j++)
    {
        Flag[0][j] = 2 * Problem[1][j];
        Flag[imax+1][j] = 1 * Problem[imax][j];
    }
    
    return 0;
}
