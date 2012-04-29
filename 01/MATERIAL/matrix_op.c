#include "matrix_op.h"
#include "helper.h"


double** add_scalar(double** M, int nrl, int nrh, int ncl, int nch, double s)
{
    int i = 0;
    int j = 0;
    double** result = matrix(nrl, nrh, ncl, nch);
    
    for(i = nrl; i <= nrh; ++i)
    {
        for(j = ncl; j <= nch; ++j)
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
    double** result = matrix(nrl, nrh, ncl, nch);
    
    for(i = nrl; i <= nrh; ++i)
    {
        for(j = ncl; j <= nch; ++j)
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
    double** result = matrix(nrl, nrh, ncl, nch);
    
    for(i = nrl; i <= nrh; ++i)
    {
        for(j = ncl; j <= nch; ++j)
        {
            result[i][j] = M[i][j] + N[i][j];
        }
    }
       
    return result;
}


double** sub_mat(double** M, double** N, int nrl, int nrh, int ncl, int nch)
{
    int i = 0;
    int j = 0;
    double** result = matrix(nrl, nrh, ncl, nch);
    
    for(i = nrl; i <= nrh; ++i)
    {
        for(j = ncl; j <= nch; ++j)
        {
            result[i][j] = M[i][j] - N[i][j];
        }
    }
       
    return result;
}
