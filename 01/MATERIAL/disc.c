#include "disc.h"
#include "helper.h"

/******** DERIVATIVES FOR EQUATION (1) START ********/

double** du2dx(double** U, double dx, int nrl, int nrh, int ncl, int nch)
{
    int i = 0;
    int j = 0;

    double firstOperand = 0;
    double secondOperand = 0;
    double alpha = 0;

    double** result = matrix(nrl, nrh, ncl, nch);
    init_matrix(result, nrl, nrh, ncl, nch, 0);

    for(i = nrl + 1; i < nrh; ++i)
    {
        for(j = ncl + 1; j < nch; ++j)
        {
            firstOperand = ( (1-alpha) / dx) *
                ( pow((U[i][j] + U[i+1][j])/2, 2) -  pow((U[i-1][j] + U[i][j])/2, 2));

            secondOperand = ( alpha / dx ) *
                ( ( (abs(U[i][j] + U[i+1][j])/2) * ((U[i][j] - U[i+1][j])/2) )  
                - ( (abs(U[i-1][j] + U[i][j])/2) * ((U[i-1][j] - U[i][j])/2) ) );

            result[i][j] = firstOperand + secondOperand;

        }
    }

   return result; 
}

double** duvdy(double** U, double **V, double dy, int nrl, int nrh, int ncl, int nch)
{
    int i = 0;
    int j = 0;

    double firstOperand = 0;
    double secondOperand = 0;
    double alpha = 0;

    double** result = matrix(nrl, nrh, ncl, nch);
    init_matrix(result, nrl, nrh, ncl, nch, 0);

    for(i = nrl + 1; i < nrh; ++i)
    {
        for(j = ncl + 1; j < nch; ++j)
        {
            firstOperand = ( ( 1- alpha ) / dy )*
            ( ( (V[i][j] + V[i+1][j]) / 2 ) * ( (U[i][j] + U[i][j+1]) / 2 )
            - ( (V[i][j-1] + V[i+1][j-1]) / 2 ) * ( (U[i][j-1] + U[i][j]) / 2 ) );

            secondOperand = ( alpha / dy )*
            ( ( (abs(V[i][j] + V[i+1][j])/2) * ((U[i][j] - U[i][j+1])/2) )
            - ( (abs(V[i][j-1] + V[i+1][j-1])/2) * ((U[i][j-1] - U[i][j])/2) ));
            result[i][j] = firstOperand + secondOperand;
        }
    }

    return result;
}


double** d2udx2(double** U, double dx, int nrl, int nrh, int ncl, int nch)
{
    int i = 0;
    int j = 0;

    double** result = matrix(nrl, nrh, ncl, nch);
    init_matrix(result, nrl, nrh, ncl, nch, 0);

    for( i = nrl + 1; i < nrh; ++i)
    {
        for( j = ncl + 1; j < nch; ++j)
        {
            result[i][j] = (U[i+1][j] - 2*U[i][j] + U[i-1][j]) / pow(dx,2);
        }
     }

     return result;
}

double** d2udy2(double** U, double dy, int nrl, int nrh, int ncl, int nch)
{
    int i = 0;
    int j = 0;

    double** result = matrix(nrl, nrh, ncl, nch);
    init_matrix(result, nrl, nrh, ncl, nch, 0);

    for( i = nrl + 1; i < nrh; ++i)
    {
        for( j = ncl + 1; j < nch; ++j)
        {
            result[i][j] = (U[i][j+1] - 2*U[i][j] + U[i][j-1]) / pow(dy,2);
        }
     }

     return result;
}

double** dpdx(double** P, double dx, int nrl, int nrh, int ncl, int nch)
{
    int i = 0;
    int j = 0;

    double** result = matrix(nrl, nrh, ncl, nch);
    init_matrix(result, nrl, nrh, ncl, nch, 0);

    for( i = nrl + 1; i < nrh; ++i)
    {
        for( j = ncl + 1; j < nch; ++j)
        {
            result[i][j]=(P[i+1][j] - P[i][j]) / dx;
        }
    }

    return result;
}

/******** DERIVATIVES FOR EQUATION (1) END ********/



/******** DERIVATIVES FOR EQUATION (2) START ********/


double** dv2dy(double** V, double dy, int nrl, int nrh, int ncl, int nch)
{
    int i = 0;
    int j = 0;

    double firstOperand = 0;
    double secondOperand = 0;
    double alpha = 0;

    double** result = matrix(nrl, nrh, ncl, nch);
    init_matrix(result, nrl, nrh, ncl, nch, 0);

    for(i = nrl + 1; i < nrh; ++i)
    {
        for(j = ncl + 1; j < nch; ++j)
        {
            firstOperand = ( (1-alpha) / dy) *
                ( pow((V[i][j] + V[i][j+1])/2, 2) -  pow((V[i][j-1] + V[i][j])/2, 2));

            secondOperand = ( alpha / dy ) *
                ( ( (abs(V[i][j] + V[i][j+1])/2) * ((V[i][j] - V[i][j+1])/2) )  
                - ( (abs(V[i][j-1] + V[i][j])/2) * ((V[i][j-1] - V[i][j])/2) ) );

            result[i][j] = firstOperand + secondOperand;

        }
    }

   return result; 
}

double** duvdx(double** U, double **V, double dx, int nrl, int nrh, int ncl, int nch)
{
    int i = 0;
    int j = 0;

    double firstOperand = 0;
    double secondOperand = 0;
    double alpha = 0;

    double** result = matrix(nrl, nrh, ncl, nch);
    init_matrix(result, nrl, nrh, ncl, nch, 0);

    for(i = nrl + 1; i < nrh; ++i)
    {
        for(j = ncl + 1; j < nch; ++j)
        {
            firstOperand = ( ( 1- alpha ) / dx )*
            ( ( (U[i][j] + U[i][j+1]) / 2 ) * ( (V[i][j] + V[i+1][j]) / 2 )
            - ( (U[i-1][j] + U[i-1][j+1]) / 2 ) * ( (V[i-1][j] + V[i][j]) / 2 ) );

            secondOperand = ( alpha / dx )*
            ( ( (abs(U[i][j] + U[i][j+1])/2) * ((V[i][j] - V[i+1][j])/2) )
            - ( (abs(U[i-1][j] + U[i-1][j+1])/2) * ((V[i-1][j] - V[i][j])/2) ));
            result[i][j] = firstOperand + secondOperand;
        }
    }

    return result;
}

double** d2vdx2(double** V, double dx, int nrl, int nrh, int ncl, int nch)
{
    int i = 0;
    int j = 0;

    double** result = matrix(nrl, nrh, ncl, nch);
    init_matrix(result, nrl, nrh, ncl, nch, 0);

    for( i = nrl + 1; i < nrh; ++i)
    {
        for( j = ncl + 1; j < nch; ++j)
        {
            result[i][j] = (V[i+1][j] - 2*V[i][j] + V[i-1][j]) / pow(dx,2);
        }
     }

     return result;
}

double** d2vdy2(double** V, double dy, int nrl, int nrh, int ncl, int nch)
{
    int i = 0;
    int j = 0;

    double** result = matrix(nrl, nrh, ncl, nch);
    init_matrix(result, nrl, nrh, ncl, nch, 0);

    for( i = nrl + 1; i < nrh; ++i)
    {
        for( j = ncl + 1; j < nch; ++j)
        {
            result[i][j] = (V[i][j+1] - 2*V[i][j] + V[i][j-1]) / pow(dy,2);
        }
     }

     return result;
}

double** dpdy(double** P, double dy, int nrl, int nrh, int ncl, int nch)
{
    int i = 0;
    int j = 0;

    double** result = matrix(nrl, nrh, ncl, nch);
    init_matrix(result, nrl, nrh, ncl, nch, 0);

    for( i = nrl + 1; i < nrh; ++i)
    {
        for( j = ncl + 1; j < nch; ++j)
        {
            result[i][j]=(P[i][j+1] - P[i][j]) / dy;
        }
    }

    return result;
}

/******** DERIVATIVES FOR EQUATION (2) END ********/













