#include "disc.h"
#include "helper.h"

double** du2dx(double** U, double dx, int nrl, int nrh, int ncl, int nch)
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
            firstOperand = ( (1-alpha) / dx) *
                ( pow((U[i][j] + U[i+1][j])/2, 2) -  pow((U[i][j] + U[i-1][j])/2, 2));

            secondOperand = ( alpha / dx ) *
                ( ( (abs(U[i][j] + U[i+1][j])/2) * ((U[i][j] - U[i+1][j])/2) )  
                - ( (abs(U[i][j] + U[i-1][j])/2) * ((U[i-1][j] - U[i][j])/2) ) );

            result[i][j] = firstOperand + secondOperand;

        }
    }

   return result; 

}

double** add_scalar(double** M, int nrl, int nrh, int ncl, int nch, double s)
{
    int i = 0;
    int j = 0;
    double** result = NULL;
    
    for(i = nrl; i < nrh; ++i)
    {
        for(j = ncl; j < nch; ++j)
        {
            result[i][j] += s;
        }
    }
    
    return result;
}

double** mult_scalar(double** M, int nrl, int nrh, int ncl, int nch, double s)
{
    int i = 0;
    int j = 0;
    double** result = NULL;
    
    for(i = nrl; i < nrh; ++i)
    {
        for(j = ncl; j < nch; ++j)
        {
            result[i][j] *= s;
        }
    }
    
    return result;
}


double** add_mat(double** M, double** N, int nrl, int nrh, int ncl, int nch)
{
    int i = 0;
    int j = 0;
    double** result = NULL;
    
    for(i = nrl; i < nrh; ++i)
    {
        for(j = ncl; j < nch; ++j)
        {
            result[i][j] = M[i][j] + M[i][j];
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

    double** result = NULL;

    for(i = nrl; i< nrh; ++i)
    {
        for(j = ncl; j < nch; ++j)
        {
            firstOperand = ( ( 1- alpha ) / dy )*
            ( ( V[i][j] + V[i+1][j] ) / 2 * ( U[i][j] + U[i][j+1] ) / 2
            - ( V[i][j-1] + V[i+1][j-1] ) / 2 * ( U[i][j] + U[i][j-1] ) / 2);

            secondOperand = ( alpha / dy )*
            ( ( (abs(V[i][j] + V[i+1][j])/2) * ((U[i][j] - U[i][j+1])/2) )
            -( (abs(V[i][j-1] + V[i+1][j-1])/2) * ((U[i][j-1] - U[i][j])/2) ));
            result[i][j] = firstOperand + secondOperand;
        }
    }

    return result;
}


double** d2udx2(double** U, double dx, int nrl, int nrh, int ncl, int nch)
{
    int i = 0;
    int j = 0;

    double** result = NULL;

    for( i = nrl; i < nrh; ++i)
    {
        for( j = ncl; j < nch; ++j)
        {
            result[i][j] = (U[i-1][j]-2*U[i][j]+U[i+1][j]) / pow(dx,2);
        }
     }


     return result;

}



double** d2udy2(double** U, double dy, int nrl, int nrh, int ncl, int nch)
{
    int i = 0;
    int j = 0;

    double** result = NULL;

    for( i = nrl; i < nrh; ++i)
    {
        for( j = ncl; j < nch; ++j)
        {
            result[i][j] = ( U[i-1][j]-2*U[i][j]+U[i+1][j]) / pow( dy,2 );
        }
    }

    return result;

}

double** dpdx(double** P, double dx, int nrl, int nrh, int ncl, int nch)
{
    int i = 0;
    int j = 0;

    double** result = NULL;

    for( i = nrl; i < nrh; ++i)
    {
        for( j = ncl; j < nch; ++j)
        {
            result[i][j]=(P[i+1][j]-P[i][j])/dx;
        }
    }

    return result;

}

















