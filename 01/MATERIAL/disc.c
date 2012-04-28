#include "disc.h"
#include "helper.h"

/******** PRIVATE DERIVATIVE FUNCTIONS START ********/
/** TODO: all the methods actually do not initialise the result before using it!! **/

double** dMds(double** M, double ds, int nrl, int nrh, int ncl, int nch)
{
    int i = 0;
    int j = 0;

    double** result = NULL;

    for( i = nrl; i < nrh; ++i)
    {
        for( j = ncl; j < nch; ++j)
        {
            result[i][j]=(M[i+1][j]-M[i][j]) / ds;
        }
    }

    return result;

}

double** dM2ds(double** M, double ds, int nrl, int nrh, int ncl, int nch)
{
    int i = 0;
    int j = 0;

    double firstOperand = 0;
    double secondOperand = 0;
    double alpha = 0;

    double** result = NULL;

    for(i = nrl; i < nrh; ++i)
    {
        for(j = ncl; j < nch; ++j)
        {
            firstOperand = ( (1-alpha) / ds) *
                ( pow((M[i][j] + M[i+1][j])/2, 2) -  pow((M[i][j] + M[i-1][j])/2, 2));

            secondOperand = ( alpha / ds ) *
                ( ( (abs(M[i][j] + M[i+1][j])/2) * ((M[i][j] - M[i+1][j])/2) )  
                - ( (abs(M[i][j] + M[i-1][j])/2) * ((M[i-1][j] - M[i][j])/2) ) );

            result[i][j] = firstOperand + secondOperand;

        }
    }

   return result; 
}


double** dMNds(double** M, double **N, double ds, int nrl, int nrh, int ncl, int nch)
{
    int i = 0;
    int j = 0;

    double firstOperand = 0;
    double secondOperand = 0;
    double alpha = 0;

    double** result = NULL;

    for(i = nrl; i< nrh; ++i)
    {
        for(j = ncl; j < nch; ++j)
        {
            firstOperand = ( ( 1- alpha ) / ds )*
            ( ( N[i][j] + N[i+1][j] ) / 2 * ( M[i][j] + M[i][j+1] ) / 2
            - ( N[i][j-1] + N[i+1][j-1] ) / 2 * ( M[i][j] + M[i][j-1] ) / 2);

            secondOperand = ( alpha / ds )*
            ( ( (abs(N[i][j] + N[i+1][j])/2) * ((M[i][j] - M[i][j+1])/2) )
            - ( (abs(N[i][j-1] + N[i+1][j-1])/2) * ((M[i][j-1] - M[i][j])/2) ));
            result[i][j] = firstOperand + secondOperand;
        }
    }

    return result;
}


double** d2Mds2(double** M, double ds, int nrl, int nrh, int ncl, int nch)
{
    int i = 0;
    int j = 0;

    double** result = NULL;

    for( i = nrl; i < nrh; ++i)
    {
        for( j = ncl; j < nch; ++j)
        {
            result[i][j] = (M[i-1][j] - 2*M[i][j] + M[i+1][j]) / pow(ds,2);
        }
     }

     return result;
}

/******** PRIVATE DERIVATIVE FUNCTIONS END ********/


/******** DERIVATIVES FOR EQUATION (1) START ********/

double** du2dx(double** U, double dx, int nrl, int nrh, int ncl, int nch)
{
   return dM2ds(U, dx, nrl, nrh, ncl, nch); 
}

double** duvdy(double** U, double **V, double dy, int nrl, int nrh, int ncl, int nch)
{
    return dMNds(U, V, dy, nrl, nrh, ncl, nch);
}


double** d2udx2(double** U, double dx, int nrl, int nrh, int ncl, int nch)
{
     return d2Mds2(U, dx, nrl, nrh, ncl, nch);
}

double** d2udy2(double** U, double dy, int nrl, int nrh, int ncl, int nch)
{
    return d2Mds2(U, dy, nrl, nrh, ncl, nch);
}

double** dpdx(double** P, double dx, int nrl, int nrh, int ncl, int nch)
{
    return dMds(P, dx, nrl, nrh, ncl, nch);
}

/******** DERIVATIVES FOR EQUATION (1) END ********/



/******** DERIVATIVES FOR EQUATION (2) START ********/


double** dv2dy(double** V, double dy, int nrl, int nrh, int ncl, int nch)
{
    return dM2ds(V, dy, nrl, nrh, ncl, nch);
}

double** duvdx(double** U, double **V, double dx, int nrl, int nrh, int ncl, int nch)
{
    return dMNds(U, V, dx, nrl, nrh, ncl, nch);
}

double** d2vdx2(double** V, double dx, int nrl, int nrh, int ncl, int nch)
{
    return d2Mds2(V, dx, nrl, nrh, ncl, nch);
}

double** d2vdy2(double** V, double dy, int nrl, int nrh, int ncl, int nch)
{
    return d2Mds2(V, dy, nrl, nrh, ncl, nch);
}

double** dpdy(double** P, double dy, int nrl, int nrh, int ncl, int nch)
{
    return dMds(P, dy, nrl, nrh, ncl, nch);
}

/******** DERIVATIVES FOR EQUATION (2) END ********/













